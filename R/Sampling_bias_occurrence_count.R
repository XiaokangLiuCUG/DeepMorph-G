# Created by Xiaokang Liu on 23/09/2026 => xkliu@cug.edu.cn
library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(icosa)
library(sp)
library(grid)

occ_path <- "D../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"
out_dir <- "./ModernGastropods/LDG"

if (!file.exists(occ_path)) stop("Occurrence file not found: ", occ_path)
if (!file.exists(morph_path)) stop("PCA_data.rds not found: ", morph_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]

pca_obj <- readRDS(morph_path)
if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else {
  morph <- as.matrix(pca_obj)
}
if (is.null(rownames(morph))) stop("Morph matrix rownames are NULL.")

occ <- occ[species %in% rownames(morph)]
if (nrow(occ) == 0) stop("No records match morphology.")

hex_num <- 5
cat(sprintf("Creating icosa grid with spacing = %d ...\n", hex_num))
hex <- hexagrid(spacing = hex_num, sf = TRUE)

coords <- cbind(occ$longitude, occ$latitude)
occ$cell_char <- locate(hex, coords)
occ <- occ[!is.na(occ$cell_char), ]
occ$cell_id <- as.integer(gsub("F", "", as.character(occ$cell_char)))

grid_summary <- occ[, .(
  occurrence_count = .N,
  species_list = list(unique(species))
), by = cell_id]
setnames(grid_summary, "cell_id", "grid_id")
grid_summary[, richness := lengths(species_list)]

cat("Converting hexagrid to sf...\n")
grid_sp <- SpPolygons(hex)
grid_sf <- st_as_sf(grid_sp)
st_crs(grid_sf) <- 4326

cell_ids_char <- sapply(slot(grid_sp, "polygons"), function(x) slot(x, "ID"))
grid_sf$grid_id <- as.integer(gsub("F", "", cell_ids_char))

grid_sf <- merge(grid_sf, grid_summary, by = "grid_id", all.x = FALSE)
grid_sf <- st_make_valid(grid_sf)

grid_sf <- grid_sf[!is.na(grid_sf$richness) & grid_sf$richness > 10, ]
cat("Number of retained grids (richness > 10):", nrow(grid_sf), "\n")

# ----------------------------  sum_var ------------------------------------
calc_diversity <- function(sp_vec, mat, pc_use = 1:10) {
  idx <- which(rownames(mat) %in% sp_vec)
  n_sp <- length(idx)
  if (n_sp < 2) return(list(n_species = n_sp, sum_var = NA))
  mat_use <- mat[idx, pc_use[pc_use <= ncol(mat)], drop = FALSE]
  if (ncol(mat_use) == 0) return(list(n_species = n_sp, sum_var = NA))
  sum_var <- sum(apply(mat_use, 2, var, na.rm = TRUE), na.rm = TRUE)
  return(list(n_species = n_sp, sum_var = sum_var))
}

sp_lists <- grid_sf$species_list
ncores <- max(1, detectCores() - 1)
if (Sys.info()["sysname"] == "Windows") {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("calc_diversity", "morph"))
  results <- parLapply(cl, sp_lists, function(sp) calc_diversity(sp, morph))
  stopCluster(cl)
} else {
  results <- mclapply(sp_lists, function(sp) calc_diversity(sp, morph), mc.cores = ncores)
}
grid_sf$sum_var <- sapply(results, function(x) x$sum_var)

# ----------------------------  sum_range ----------------------------------
calc_sum_range <- function(sp_vec, mat, pc_use = 1:10) {
  idx <- which(rownames(mat) %in% sp_vec)
  if (length(idx) < 2) return(NA)
  mat_use <- mat[idx, pc_use[pc_use <= ncol(mat)], drop = FALSE]
  if (ncol(mat_use) == 0) return(NA)
  range_sum <- sum(apply(mat_use, 2, function(x) diff(range(x, na.rm = TRUE))), na.rm = TRUE)
  return(range_sum)
}

scores_10 <- morph[, 1:min(10, ncol(morph)), drop = FALSE]
if (Sys.info()["sysname"] == "Windows") {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("calc_sum_range", "scores_10"))
  res <- parLapply(cl, sp_lists, function(sp) calc_sum_range(sp, scores_10))
  stopCluster(cl)
} else {
  res <- mclapply(sp_lists, function(sp) calc_sum_range(sp, scores_10), mc.cores = ncores)
}
grid_sf$sum_range <- unlist(res)

# ----------------------------  local_jaccard_mean -------------------------
all_species <- unique(unlist(grid_sf$species_list))
mat <- matrix(0, nrow = nrow(grid_sf), ncol = length(all_species))
rownames(mat) <- grid_sf$grid_id
colnames(mat) <- all_species
for (i in seq_len(nrow(grid_sf))) {
  sp <- grid_sf$species_list[[i]]
  mat[i, sp] <- 1
}

all_neighbors <- st_touches(grid_sf, grid_sf)
grid_sf$local_jaccard_mean <- NA

for (i in seq_len(nrow(grid_sf))) {
  neigh <- all_neighbors[[i]]
  if (length(neigh) == 0) next
  valid_neigh <- neigh[neigh %in% 1:nrow(grid_sf)]
  if (length(valid_neigh) == 0) next
  jaccard_vals <- sapply(valid_neigh, function(j) {
    a <- mat[i, ]
    b <- mat[j, ]
    inter <- sum(a & b)
    union <- sum(a | b)
    if (union == 0) return(NA)
    1 - inter / union
  })
  grid_sf$local_jaccard_mean[i] <- mean(jaccard_vals, na.rm = TRUE)
}

# ----------------------------  MSp ----------------------------------------
calc_MSp <- function(sp_vec, mat, pc_use = 1:10) {
  idx <- which(rownames(mat) %in% sp_vec)
  if (length(idx) < 2) return(NA)
  coords <- mat[idx, pc_use[pc_use <= ncol(mat)], drop = FALSE]
  if (ncol(coords) == 0) return(NA)
  centroid <- colMeans(coords)
  dist_cent <- sqrt(rowSums((coords - centroid)^2))
  return(mean(dist_cent))
}

if (Sys.info()["sysname"] == "Windows") {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("calc_MSp", "morph"))
  res_ms <- parLapply(cl, sp_lists, function(sp) calc_MSp(sp, morph))
  stopCluster(cl)
} else {
  res_ms <- mclapply(sp_lists, function(sp) calc_MSp(sp, morph), mc.cores = ncores)
}
grid_sf$MSp <- unlist(res_ms)

metrics <- c("richness", "sum_var", "sum_range", "local_jaccard_mean", "MSp")

metric_labels <- c(
  "richness"           = "Species richness",
  "sum_var"            = "Morph. disparity (sum var)",
  "sum_range"          = "Morph. disparity (sum range)",
  "local_jaccard_mean" = "Local beta diversity",
  "MSp"                = "Morph. specialization"
)

log_y_metrics <- c("richness")   


dat_raw <- grid_sf %>%
  st_set_geometry(NULL) %>%
  select(occurrence_count, all_of(metrics))


dat_long <- dat_raw %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "Metric_raw", values_to = "value") %>%
  filter(!is.na(occurrence_count), !is.na(value), occurrence_count > 0) %>%
  mutate(Metric = factor(Metric_raw,
                         levels = names(metric_labels),
                         labels = metric_labels))

# -----------------------------------------------------------------------------
# Per-metric: fit + R² + p-value + fitted curve
# -----------------------------------------------------------------------------
fit_metric <- function(m_raw, df_all) {
  df <- df_all %>% filter(Metric_raw == m_raw)
  if (nrow(df) < 3) return(NULL)
  
  use_log_y <- m_raw %in% log_y_metrics
  
  if (use_log_y) {
    fit <- lm(log10(value) ~ log10(occurrence_count), data = df)
  } else {
    fit <- lm(value ~ log10(occurrence_count), data = df)
  }
  s <- summary(fit)
  
  xseq <- seq(min(df$occurrence_count), max(df$occurrence_count),
              length.out = 200)
  if (use_log_y) {
    yhat <- 10^predict(fit, newdata = data.frame(occurrence_count = xseq))
  } else {
    yhat <- predict(fit, newdata = data.frame(occurrence_count = xseq))
  }
  
  p_val   <- s$coefficients[2, "Pr(>|t|)"]
  p_label <- ifelse(p_val < 2.2e-16, "< 2.2e-16", sprintf("%.2g", p_val))
  
  list(
    stat = tibble(
      Metric_raw = m_raw,
      Metric     = metric_labels[[m_raw]],
      n          = nrow(df),
      R2         = s$r.squared,
      adj_R2     = s$adj.r.squared,
      p          = p_val,
      slope      = s$coefficients[2, "Estimate"],
      p_label    = p_label,
      label      = sprintf("R\u00B2 = %.3f\np %s\nn = %d",
                           s$r.squared, p_label, nrow(df))
    ),
    line = tibble(occurrence_count = xseq, yhat = yhat),
    data = df,
    log_y = use_log_y
  )
}

fits <- lapply(metrics, fit_metric, df_all = dat_long)
names(fits) <- metrics

fit_res <- bind_rows(lapply(fits, "[[", "stat"))
print(fit_res)
fwrite(as.data.frame(fit_res),
       file.path(out_dir, "OccurrenceCount_vs_metrics_R2.csv"))


make_metric_plot <- function(fit_obj, show_xlab = TRUE, show_ylab = TRUE) {
  if (is.null(fit_obj)) return(NULL)
  st  <- fit_obj$stat
  ln  <- fit_obj$line
  df  <- fit_obj$data
  use_log_y <- fit_obj$log_y
  
  p <- ggplot(df, aes(x = occurrence_count, y = value)) +
    geom_point(alpha = 0.35, size = 0.7, color = "steelblue") +
    geom_line(data = ln, aes(x = occurrence_count, y = yhat),
              color = "firebrick", linewidth = 0.9) +
    annotate("text",
             x = Inf, y = Inf,
             label = st$label,
             hjust = 1.05, vjust = 1.1,
             size = 2.8, color = "black") +
    scale_x_log10(
      breaks = scales::log_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +
    ggtitle(st$Metric) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.margin = margin(4, 6, 4, 4)
    )
  
  if (use_log_y) {
    p <- p + scale_y_log10(
      breaks = scales::log_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    )
  }
  
  p <- p + labs(
    x = if (show_xlab) "Occurrence count (log10 scale)" else NULL,
    y = if (show_ylab) "Metric value" else NULL
  )
  p
}

p_richness <- make_metric_plot(fits[["richness"]])
p_sumvar   <- make_metric_plot(fits[["sum_var"]])
p_sumrange <- make_metric_plot(fits[["sum_range"]])
p_jaccard  <- make_metric_plot(fits[["local_jaccard_mean"]])
p_msp      <- make_metric_plot(fits[["MSp"]])



pdf(file.path(out_dir, "OccurrenceCount_vs_metrics_R2.pdf"),
    width = 11, height = 13)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 3, ncol = 2,
  heights = unit(c(1, 1, 1), "null"),
  widths  = unit(c(1, 1), "null")
)))


grid.text(
  "Relationship between occurrence count and five diversity metrics",
  y = unit(0.99, "npc"),
  gp = gpar(fontface = "bold", fontsize = 13)
)
grid.text(
  "Raw grid-cell values (no resampling); richness on log-log scale, others on log10(x); richness y-axis on log10 scale",
  y = unit(0.965, "npc"),
  gp = gpar(fontsize = 8, col = "grey30")
)

# Panel 1: richness (spanning two columns)
print(p_richness,
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# Panels 2 and 3
print(p_sumvar,
      vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_sumrange,
      vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# Panels 4 and 5
print(p_jaccard,
      vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p_msp,
      vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

popViewport()
dev.off()

cat("Done. Output:\n",
    file.path(out_dir, "OccurrenceCount_vs_metrics_R2.pdf"), "\n",
    file.path(out_dir, "OccurrenceCount_vs_metrics_R2.csv"), "\n")


vars_export <- c("grid_id", "occurrence_count", "richness",
                 "sum_var", "sum_range", "local_jaccard_mean", "MSp")

grid_values <- grid_sf %>%
  st_set_geometry(NULL) %>%
  select(all_of(vars_export)) %>%
  arrange(grid_id)

# Optional: also save longitude/latitude for downstream spatial matching
grid_coords <- grid_sf %>%
  mutate(
    lon = st_coordinates(st_point_on_surface(geometry))[, "X"],
    lat = st_coordinates(st_point_on_surface(geometry))[, "Y"]
  ) %>%
  st_set_geometry(NULL) %>%
  select(grid_id, lon, lat)

grid_values <- grid_values %>%
  left_join(grid_coords, by = "grid_id") %>%
  select(grid_id, lon, lat,
         occurrence_count, richness,
         sum_var, sum_range, local_jaccard_mean, MSp)

fwrite(as.data.frame(grid_values),
       file.path(out_dir, "Gridcell_six_variables.csv"))

cat("Saved:", file.path(out_dir, "Gridcell_six_variables.csv"),
    "  nrow =", nrow(grid_values), "\n")

# Quick view of the first few rows
print(head(grid_values))