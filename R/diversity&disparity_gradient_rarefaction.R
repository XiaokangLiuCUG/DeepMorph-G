# Created by Xiaokang Liu on 20/10/2026 => xkliu@cug.edu.cn
# Direct partitioning by 2.5° latitude bands + rarefaction sampling
# Metrics: richness, sum_var, sum_range, local_jaccard_mean, MSp

library(data.table)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# ---------------------------- Paths -------------------------------------------
occ_path   <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"
out_dir    <- "./ModernGastropods"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Read data ---------------------------------------
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]

pca_obj <- readRDS(morph_path)
morph <- if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) pca_obj$x else as.matrix(pca_obj)
if (is.null(rownames(morph))) stop("Morph matrix rownames are NULL.")

occ <- occ[species %in% rownames(morph)]
if (nrow(occ) == 0) stop("No records match morphology.")

# ---------------------------- Partition by 2.5° latitude bands (no hexagonal grid) -----------
lat_breaks <- seq(-90, 90, by = 2.5)
occ[, lat_bin := cut(latitude, breaks = lat_breaks, include.lowest = TRUE)]
occ <- occ[!is.na(lat_bin)]

lat_mids <- lat_breaks[-1] - diff(lat_breaks) / 2
names(lat_mids) <- levels(occ$lat_bin)

# ---------------------------- Count occurrences per band -------------------------
occ_per_band <- occ[, .(n_occ = .N, n_species = uniqueN(species)), by = lat_bin]
occ_per_band[, lat_mid := lat_mids[as.character(lat_bin)]]
setorder(occ_per_band, lat_mid)

cat("\n=== Occurrence count per 2.5° latitude band ===\n")
print(occ_per_band, n = Inf)

cat("\n=== Quantiles of occurrence counts ===\n")
print(quantile(occ_per_band$n_occ, c(0, 0.1, 0.2, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE))

fwrite(occ_per_band, file.path(out_dir, "Occurrence_count_by_band.csv"))

# ---------------------------- Choose rarefaction sample size N_rare ----------------------------############################################################
q20 <- as.numeric(quantile(occ_per_band$n_occ, 0.2, na.rm = TRUE))
N_rare <- max(500, floor(q20))

n_bands_total  <- nrow(occ_per_band)
n_bands_usable <- sum(occ_per_band$n_occ >= N_rare)
cat(sprintf("\nN_rare = %d (20th percentile = %.0f)\n", N_rare, q20))
cat(sprintf("Bands eligible for rarefaction: %d / %d (%.1f%%)\n",
            n_bands_usable, n_bands_total, 100 * n_bands_usable / n_bands_total))

if (n_bands_usable / n_bands_total < 0.8) {
  q10 <- as.numeric(quantile(occ_per_band$n_occ, 0.1, na.rm = TRUE))
  N_rare <- max(50, floor(q10))
  n_bands_usable <- sum(occ_per_band$n_occ >= N_rare)
  cat(sprintf("Usable rate < 80%%, automatically lowering to 10th percentile: N_rare = %d, usable %d / %d\n",
              N_rare, n_bands_usable, n_bands_total))
}

#

calc_metrics_from_occ <- function(dt_sub, morph, pc_use = 1:10) {
  sp <- unique(dt_sub$species)
  n_sp <- length(sp)
  idx <- which(rownames(morph) %in% sp)
  if (length(idx) < 2) {
    return(list(richness = n_sp, sum_var = NA, sum_range = NA,
                local_jaccard_mean = NA, MSp = NA))
  }
  mat_use <- morph[idx, pc_use[pc_use <= ncol(morph)], drop = FALSE]
  sum_var   <- sum(apply(mat_use, 2, var, na.rm = TRUE), na.rm = TRUE)
  sum_range <- sum(apply(mat_use, 2, function(x) diff(range(x, na.rm = TRUE))), na.rm = TRUE)
  centroid  <- colMeans(mat_use)
  MSp       <- mean(sqrt(rowSums((mat_use - centroid)^2)))
  list(richness = n_sp, sum_var = sum_var, sum_range = sum_range,
       local_jaccard_mean = NA, MSp = MSp)
}

# ---------------------------- Single-band rarefaction function -------------------------------
rarefy_one_band <- function(b, dt, N_rare, R, morph) {
  dt_b <- dt[lat_bin == b]
  if (nrow(dt_b) < N_rare) return(NULL)
  res <- replicate(R, {
    idx <- sample.int(nrow(dt_b), N_rare, replace = FALSE)
    calc_metrics_from_occ(dt_b[idx], morph)
  }, simplify = FALSE)
  df <- bind_rows(lapply(res, as.data.frame))
  df$lat_bin <- b
  df
}

# ---------------------------- Rarefaction + R replicates
R <- 1000
set.seed(42)   # reproducible

usable_bands <- occ_per_band[n_occ >= N_rare, lat_bin]

rarefied_list <- lapply(usable_bands, function(b) {
  rarefy_one_band(b, occ, N_rare, R, morph)
})
rarefied_df <- bind_rows(rarefied_list)

# Remove bands with fewer than 2 unique species that cause all-NA metrics
rarefied_df <- rarefied_df %>%
  filter(!is.na(richness))

cat(sprintf("\nRarefaction complete: %d bands, %d replicates per band\n",
            length(unique(rarefied_df$lat_bin)), R))

# ---------------------------- Summarize
rarefied_summary <- rarefied_df %>%
  group_by(lat_bin) %>%
  summarise(
    n_rep = n(),
    richness_mean  = mean(richness,  na.rm = TRUE),
    richness_lower = quantile(richness,  0.025, na.rm = TRUE),
    richness_upper = quantile(richness,  0.975, na.rm = TRUE),
    sum_var_mean   = mean(sum_var,   na.rm = TRUE),
    sum_var_lower  = quantile(sum_var,   0.025, na.rm = TRUE),
    sum_var_upper  = quantile(sum_var,   0.975, na.rm = TRUE),
    sum_range_mean = mean(sum_range, na.rm = TRUE),
    sum_range_lower= quantile(sum_range, 0.025, na.rm = TRUE),
    sum_range_upper= quantile(sum_range, 0.975, na.rm = TRUE),
    local_jaccard_mean  = mean(local_jaccard_mean,  na.rm = TRUE),
    local_jaccard_lower = quantile(local_jaccard_mean, 0.025, na.rm = TRUE),
    local_jaccard_upper = quantile(local_jaccard_mean, 0.975, na.rm = TRUE),
    MSp_mean  = mean(MSp,  na.rm = TRUE),
    MSp_lower = quantile(MSp,  0.025, na.rm = TRUE),
    MSp_upper = quantile(MSp,  0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lat_mid = lat_mids[as.character(lat_bin)]) %>%
  filter(!is.na(lat_mid)) %>%
  arrange(lat_mid)

# Convert NaN (all-NA metrics) to NA for plotting
rarefied_summary <- rarefied_summary %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

fwrite(as.data.table(rarefied_summary),
       file.path(out_dir, "Rarefied_metrics_summary_by_band.csv"))

# ---------------------------- Convert to long format + plotting -------------------------------
metric_labels <- c(
  "richness"      = "Species richness",
  "sum_var"       = "Morph. disparity (sum var)",
  "sum_range"     = "Morph. disparity (sum range)",
  "local_jaccard" = "Local beta diversity",
  "MSp"           = "Morph. specialization"
)

# Detect metrics that are all NA (local_jaccard is always NA in no-grid mode)
available_metrics <- names(metric_labels)
if (all(is.na(rarefied_summary$local_jaccard_mean))) {
  available_metrics <- setdiff(available_metrics, "local_jaccard")
  cat("\n[Note] local_jaccard_mean cannot be computed in no-grid mode and has been removed from the plot.\n")
}

lat_long <- rarefied_summary %>%
  select(lat_bin, lat_mid,
         richness_mean, richness_lower, richness_upper,
         sum_var_mean, sum_var_lower, sum_var_upper,
         sum_range_mean, sum_range_lower, sum_range_upper,
         local_jaccard_mean, local_jaccard_lower, local_jaccard_upper,
         MSp_mean, MSp_lower, MSp_upper) %>%
  pivot_longer(
    cols = -c(lat_bin, lat_mid),
    names_to = c("Metric", ".value"),
    names_pattern = "(.*)_(mean|lower|upper)"
  ) %>%
  filter(Metric %in% available_metrics) %>%
  mutate(Metric = factor(Metric,
                         levels = available_metrics,
                         labels = metric_labels[available_metrics]))

p <- ggplot(lat_long, aes(x = lat_mid, y = mean, color = Metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Metric),
              alpha = 0.25, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 0.8, na.rm = TRUE) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    x = "Latitude (°)",
    y = "Mean value",
    subtitle = sprintf("Rarefied to N = %d occurrences per band, R = %d replicates, 95%% CI",
                       N_rare, R)
  ) +
  ggtitle("Latitudinal gradients of five metrics (rarefied by sampling effort)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()
  ) +
  scale_x_continuous(breaks = seq(-90, 90, by = 15),
                     limits = c(-90, 90), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-90, 90))

ggsave(file.path(out_dir, "Latitudinal_Gradients_5metrics_rarefied.pdf"),
       plot = p, width = 10, height = 12, dpi = 300)
print(p)