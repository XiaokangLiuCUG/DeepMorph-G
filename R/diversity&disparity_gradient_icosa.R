# Created by Xiaokang Liu on 30/08/2026 => xkliu@cug.edu.cn
# Output：Latitudinal_Gradients_six_metrics_with_CI.pdf
# Metrics：richness, occurrence_count, sum_var, sum_range, local_jaccard_mean, MSp
# =============================================================================

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


occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"
out_dir <- "./ModernGastropods/LGD"

if (!file.exists(occ_path)) stop("Occurrence file not found: ", occ_path)
if (!file.exists(morph_path)) stop("PCA_data.rds not found: ", morph_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- data ----------------------------------------
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

# ---------------------------- cell（icosa）-------------------
hex_num <- 5   # （ 5, 6, 10）
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

# ======================== richness > 10 =================
grid_sf <- grid_sf[!is.na(grid_sf$richness) & grid_sf$richness > 10, ]
cat("保留网格数（richness > 10）:", nrow(grid_sf), "\n")


grid_sf <- grid_sf %>% 
  mutate(lat = st_coordinates(st_point_on_surface(geometry))[, "Y"])

# ---------------------------- sum_var-----------------------------
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

# ---------------------------- sum_range------------------
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

# ---------------------------- local_jaccard_mean -------------------------
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

# ---------------------------- MSp（形态特化程度）-------------------------
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

# ---------------------------- （2.5°） + Bootstrap 置信区间 -----------
lat_breaks <- seq(-90, 90, by = 2.5)
grid_sf$lat_bin <- cut(grid_sf$lat, breaks = lat_breaks, include.lowest = TRUE)
lat_mids <- lat_breaks[-1] - diff(lat_breaks)/2
names(lat_mids) <- levels(grid_sf$lat_bin)

B <- 1000
boot_ci <- function(x, B = 1000, conf = 0.95) {
  if (length(x) < 2) return(tibble(mean = NA, lower = NA, upper = NA))
  boot_means <- replicate(B, {
    sample(x, size = length(x), replace = TRUE) %>% mean(na.rm = TRUE)
  })
  alpha <- (1 - conf) / 2
  tibble(
    mean = mean(x, na.rm = TRUE),
    lower = quantile(boot_means, alpha, na.rm = TRUE),
    upper = quantile(boot_means, 1 - alpha, na.rm = TRUE)
  )
}

lat_ci <- grid_sf %>%
  st_set_geometry(NULL) %>%
  filter(!is.na(lat_bin)) %>%
  group_by(lat_bin) %>%
  summarise(
    richness_res = list(boot_ci(richness, B = B)),
    occ_count_res = list(boot_ci(occurrence_count, B = B)),
    sum_var_res = list(boot_ci(sum_var, B = B)),
    sum_range_res = list(boot_ci(sum_range, B = B)),
    local_jaccard_res = list(boot_ci(local_jaccard_mean, B = B)),
    MSp_res = list(boot_ci(MSp, B = B)),
    .groups = "drop"
  ) %>%
  mutate(
    richness_mean = map_dbl(richness_res, "mean"),
    richness_lower = map_dbl(richness_res, "lower"),
    richness_upper = map_dbl(richness_res, "upper"),
    occ_count_mean = map_dbl(occ_count_res, "mean"),
    occ_count_lower = map_dbl(occ_count_res, "lower"),
    occ_count_upper = map_dbl(occ_count_res, "upper"),
    sum_var_mean = map_dbl(sum_var_res, "mean"),
    sum_var_lower = map_dbl(sum_var_res, "lower"),
    sum_var_upper = map_dbl(sum_var_res, "upper"),
    sum_range_mean = map_dbl(sum_range_res, "mean"),
    sum_range_lower = map_dbl(sum_range_res, "lower"),
    sum_range_upper = map_dbl(sum_range_res, "upper"),
    local_jaccard_mean = map_dbl(local_jaccard_res, "mean"),
    local_jaccard_lower = map_dbl(local_jaccard_res, "lower"),
    local_jaccard_upper = map_dbl(local_jaccard_res, "upper"),
    MSp_mean = map_dbl(MSp_res, "mean"),
    MSp_lower = map_dbl(MSp_res, "lower"),
    MSp_upper = map_dbl(MSp_res, "upper")
  ) %>%
  mutate(lat_mid = lat_mids[as.character(lat_bin)]) %>%
  filter(!is.na(lat_mid))


metric_labels <- c(
  "richness" = "Species richness",
  "occ_count" = "Occurrence count",
  "sum_var" = "Morph. disparity (sum var)",
  "sum_range" = "Morph. disparity (sum range)",
  "local_jaccard" = "Local beta diversity",
  "MSp" = "Morph. specialization"
)

lat_ci_long <- lat_ci %>%
  select(lat_bin, lat_mid, 
         richness_mean, richness_lower, richness_upper,
         occ_count_mean, occ_count_lower, occ_count_upper,
         sum_var_mean, sum_var_lower, sum_var_upper,
         sum_range_mean, sum_range_lower, sum_range_upper,
         local_jaccard_mean, local_jaccard_lower, local_jaccard_upper,
         MSp_mean, MSp_lower, MSp_upper) %>%
  pivot_longer(
    cols = -c(lat_bin, lat_mid),
    names_to = c("Metric", ".value"),
    names_pattern = "(.*)_(mean|lower|upper)"
  ) %>%
  mutate(Metric = factor(Metric,
                         levels = names(metric_labels),
                         labels = metric_labels))

# ---------------------------- plotting----------------------------
p_lat_ci <- ggplot(lat_ci_long, aes(x = lat_mid, y = mean, color = Metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Metric), 
              alpha = 0.25, color = NA) +
  geom_line(size = 1) + 
  geom_point(size = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(x = "Latitude (°)", y = "Mean value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank()
  ) +
  ggtitle("Latitudinal gradients of six diversity metrics with 95% CI (bootstrap)") +
  scale_x_continuous(
    breaks = seq(-90, 90, by = 15),
    limits = c(-90, 90),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(-90, 90))

ggsave(file.path(out_dir, "Latitudinal_Gradients_six_metrics_with_CI_icosa.pdf"),
       plot = p_lat_ci, width = 10, height = 14, dpi = 300)

print(p_lat_ci)



# =============================================================================
# save raw data
# =============================================================================

original_metrics <- grid_sf %>%
  st_set_geometry(NULL) %>%
  select(
    grid_id,
    lat,           # 网格质心纬度
    lat_bin,       # 2.5° 纬度带
    richness,            # 物种丰富度（原始值）
    occurrence_count,    # 出现记录数（原始值）
    sum_var,             # 形态方差和（原始值）
    sum_range,           # 形态极差和（原始值）
    local_jaccard_mean,  # 局部 beta 多样性（原始值）
    MSp                  # 形态特化程度（原始值）
  ) %>%
  arrange(lat)


print(head(original_metrics))

cat("\n原始数据总行数（网格数）：", nrow(original_metrics), "\n")

# --- 3. save CSV ---
out_csv <- file.path(out_dir, "Six_metrics_original_by_grid.csv")
fwrite(original_metrics, out_csv)
cat("\n原始数据已保存至：", out_csv, "\n")


original_by_latbin <- original_metrics %>%
  filter(!is.na(lat_bin))

out_csv_latbin <- file.path(out_dir, "Six_metrics_original_by_latbin.csv")
fwrite(original_by_latbin, out_csv_latbin)
cat("按纬度带分组的原始数据已保存至：", out_csv_latbin, "\n")


original_long <- original_metrics %>%
  filter(!is.na(lat_bin)) %>%
  select(lat_bin, lat, richness, occurrence_count, sum_var,
         sum_range, local_jaccard_mean, MSp) %>%
  pivot_longer(
    cols = c(richness, occurrence_count, sum_var,
             sum_range, local_jaccard_mean, MSp),
    names_to = "Metric",
    values_to = "value"
  )

out_csv_long <- file.path(out_dir, "Six_metrics_original_long.csv")
fwrite(original_long, out_csv_long)
cat("长格式原始数据已保存至：", out_csv_long, "\n")








# 六个指标全部 log10
lat_ci_long <- lat_ci_long %>%
  mutate(
    mean_plot  = log10(pmax(mean,  1e-10)),
    lower_plot = log10(pmax(lower, 1e-10)),
    upper_plot = log10(pmax(upper, 1e-10))
  )

p_lat_ci <- ggplot(lat_ci_long, aes(x = lat_mid, y = mean_plot, color = Metric)) +
  geom_ribbon(aes(ymin = lower_plot, ymax = upper_plot, fill = Metric),
              alpha = 0.25, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_y_continuous(
    labels = scales::math_format(10^.x)
  ) +
  labs(x = "Latitude (°)", y = "Mean value (log10 scale)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank()
  ) +
  ggtitle("Latitudinal gradients of six diversity metrics with 95% CI (bootstrap)") +
  scale_x_continuous(
    breaks = seq(-90, 90, by = 15),
    limits = c(-90, 90),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(-90, 90))

ggsave(file.path(out_dir, "Latitudinal_Gradients_six_metrics_with_CI_icosa_log10.pdf"),
       plot = p_lat_ci, width = 10, height = 14, dpi = 300)

print(p_lat_ci)
