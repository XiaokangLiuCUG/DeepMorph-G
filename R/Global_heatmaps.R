# Created by Xiaokang Liu on 15/7/2026 => xkliu@cug.edu.cn
# =============================================================================
# Gastropod diversity and morphological analysis - icosa equal-area hexagon grid + latitudinal band fitting
# =============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(parallel)
library(viridis)
library(icosa)
library(sp)


occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"

if (!file.exists(occ_path)) stop("Occurrence file not found: ", occ_path)
if (!file.exists(morph_path)) stop("Morphology file not found: ", morph_path)

occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
cat("Total occurrence records:", nrow(occ), "\n")

pca_obj <- readRDS(morph_path)
if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else {
  morph <- as.matrix(pca_obj)
}
if (is.null(rownames(morph))) stop("Morph matrix rownames are NULL.")
cat("Morph matrix dimensions:", dim(morph), "\n")

occ <- occ[species %in% rownames(morph)]
cat("Records after morphology matching:", nrow(occ), "\n")
if (nrow(occ) == 0) stop("No records match morphology.")


hex_num <- 5   # changed to 5 (finer); increase slightly if errors persist
cat(sprintf("Creating hexagrid with spacing = %d ...\n", hex_num))
hex <- hexagrid(spacing = hex_num, sf = TRUE)
cat(sprintf("Total hexagon cells: %d\n", length(hex)))

occ$cell <- locate(hex, occ[, c("longitude", "latitude")])
occ <- occ[!is.na(occ$cell), ]
cat(sprintf("Records assigned to cells: %d\n", nrow(occ)))

setDT(occ)
grid_summary <- occ[, .(
  occurrence_count = .N,
  species_list = list(unique(species))
), by = cell]
grid_summary[, richness := lengths(species_list)]
grid_summary[, cell := as.character(cell)]

grid_plot <- grid_summary[richness > 10, ]
cat("Grids with richness > 10:", nrow(grid_plot), "\n")
if (nrow(grid_plot) == 0) stop("No cells with richness > 10. Try lowering hex_num or threshold.")

cat("Converting hexagrid to SpatialPolygons ...\n")
grid_sp <- SpPolygons(hex)
cell_ids <- sapply(slot(grid_sp, "polygons"), function(x) slot(x, "ID"))

grid_sf <- st_as_sf(grid_sp)
st_crs(grid_sf) <- 4326
grid_sf$cell <- cell_ids

# ★ Key fix: use st_make_valid to repair invalid geometries
grid_sf <- st_make_valid(grid_sf)

# Merge statistics
grid_sf <- merge(grid_sf, grid_plot, by = "cell", all.x = FALSE)
cat("Valid grids after filtering (richness > 10):", nrow(grid_sf), "\n")
if (nrow(grid_sf) == 0) stop("Merging failed: no matching cells. Check cell names.")

# ---------------------------- 7. Compute sum_var and sum_range --------------------
scores <- morph
n_pcs <- min(10, ncol(scores))
scores_10 <- scores[, 1:n_pcs, drop = FALSE]

calc_disparity <- function(sp_list, scores_mat) {
  if (length(sp_list) < 2) return(c(sum_var = NA, sum_range = NA))
  idx <- which(rownames(scores_mat) %in% sp_list)
  if (length(idx) < 2) return(c(sum_var = NA, sum_range = NA))
  mat <- scores_mat[idx, , drop = FALSE]
  sum_var   <- sum(apply(mat, 2, var, na.rm = TRUE), na.rm = TRUE)
  sum_range <- sum(apply(mat, 2, function(x) diff(range(x, na.rm = TRUE))), na.rm = TRUE)
  return(c(sum_var = sum_var, sum_range = sum_range))
}

grid_sf$sum_var <- NA
grid_sf$sum_range <- NA

sp_lists <- grid_sf$species_list
if (nrow(grid_sf) > 100) {
  ncores <- max(1, detectCores() - 1)
  cat("Using", ncores, "cores for parallel computation...\n")
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("calc_disparity", "scores_10"))
    results <- parLapply(cl, sp_lists, function(sp) calc_disparity(sp, scores_10))
    stopCluster(cl)
  } else {
    results <- mclapply(sp_lists, function(sp) calc_disparity(sp, scores_10), mc.cores = ncores)
  }
  for (i in 1:nrow(grid_sf)) {
    grid_sf$sum_var[i] <- results[[i]]["sum_var"]
    grid_sf$sum_range[i] <- results[[i]]["sum_range"]
  }
} else {
  for (i in 1:nrow(grid_sf)) {
    res <- calc_disparity(grid_sf$species_list[[i]], scores_10)
    grid_sf$sum_var[i] <- res["sum_var"]
    grid_sf$sum_range[i] <- res["sum_range"]
  }
}
cat("Disparity calculation completed.\n")

# ★ After repairing geometry, re-extract centroids (use st_point_on_surface to avoid invalid geometry issues)
cat("Computing centroids...\n")
grid_sf <- st_make_valid(grid_sf)   # ensure validity again
centroids <- st_point_on_surface(grid_sf$geometry)   # more robust
centroid_coords <- st_coordinates(centroids)
grid_sf$centroid_lat <- centroid_coords[, "Y"]
grid_sf$abs_centroid_lat <- abs(grid_sf$centroid_lat)

# Assign latitudinal bands
grid_sf$lat_band <- cut(grid_sf$abs_centroid_lat,
                        breaks = c(0, 30, 45, 60, 90),
                        include.lowest = TRUE,
                        labels = c("Lower", "Second", "Third", "Upper"))
band_colors <- c("Lower"   = "#000080",
                 "Second"  = "#33a02c",
                 "Third"   = "#ff7f00",
                 "Upper"   = "#6a3d9a")

# ---------------------------- 8. Plot heatmaps (Equirectangular projection) --------------
world <- ne_countries(scale = "medium", returnclass = "sf")

plot_heatmap <- function(data, var_name, title_prefix = "Global distribution of") {
  vals <- data[[var_name]]
  vals_pos <- vals[vals > 0 & !is.na(vals)]
  if (length(vals_pos) == 0) {
    warning(paste("No positive values for", var_name, "- skipping plot."))
    return(NULL)
  }
  colours <- viridis(100, option = "plasma")
  p <- ggplot() +
    geom_sf(data = world, fill = "gray50", color = "gray70", size = 0.2) +
    geom_sf(data = data,
            aes(fill = ifelse(.data[[var_name]] > 0, .data[[var_name]], NA)),
            color = NA, alpha = 1) +
    scale_fill_gradientn(
      colours = colours,
      na.value = NA,
      breaks = pretty(vals_pos, n = 5),
      labels = scales::comma,
      #guide = "none",
      name = var_name
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = paste(title_prefix, var_name))
  return(p)
}

p_occurrence <- plot_heatmap(grid_sf, "occurrence_count")
p_richness   <- plot_heatmap(grid_sf, "richness")
p_sumvar     <- plot_heatmap(grid_sf, "sum_var")
p_sumrange   <- plot_heatmap(grid_sf, "sum_range")

if (!is.null(p_occurrence)) {
  ggsave("Global_occurrence_count_heatmap_icosa0.pdf", p_occurrence, width = 12, height = 6)
  print(p_occurrence)
}
if (!is.null(p_richness)) {
  ggsave("Global_richness_heatmap_icosa0.pdf", p_richness, width = 12, height = 6)
  print(p_richness)
}
if (!is.null(p_sumvar)) {
  ggsave("Global_sum_var_heatmap_icosa0.pdf", p_sumvar, width = 12, height = 6)
  print(p_sumvar)
}
if (!is.null(p_sumrange)) {
  ggsave("Global_sum_range_heatmap_icosa0.pdf", p_sumrange, width = 12, height = 6)
  print(p_sumrange)
}

# ---------------------------- Latitudinal band fitting scatter plots (continuation) ------------------
# This section follows the code you provided: for each latitudinal band, fit log10(y) ~ log10(x)
bands <- levels(grid_sf$lat_band)

# Store label positions and data
label_data_var <- data.frame()
label_data_range <- data.frame()

for (b in bands) {
  sub <- grid_sf[grid_sf$lat_band == b, ]
  sub <- sub[!is.na(sub$richness) & !is.na(sub$sum_var) & !is.na(sub$sum_range) & 
               sub$richness > 0 & sub$sum_var > 0 & sub$sum_range > 0, ]
  if (nrow(sub) < 3) next
  
  # sum_var
  fit_var <- lm(log10(sum_var) ~ log10(richness), data = sub)
  cf_var <- coef(fit_var)
  r2_var <- summary(fit_var)$r.squared
  label_var <- sprintf("log10(y) = %.3f*log10(x) + %.3f, R²=%.3f", cf_var[2], cf_var[1], r2_var)
  
  # sum_range
  fit_range <- lm(log10(sum_range) ~ log10(richness), data = sub)
  cf_range <- coef(fit_range)
  r2_range <- summary(fit_range)$r.squared
  label_range <- sprintf("log10(y) = %.3f*log10(x) + %.3f, R²=%.3f", cf_range[2], cf_range[1], r2_range)
  
  x_pos <- max(sub$richness, na.rm = TRUE)
  y_pos_var <- max(sub$sum_var, na.rm = TRUE) * 0.9
  y_pos_range <- max(sub$sum_range, na.rm = TRUE) * 0.9
  
  label_data_var <- rbind(label_data_var, data.frame(
    band = b, x = x_pos, y = y_pos_var, label = label_var, stringsAsFactors = FALSE
  ))
  label_data_range <- rbind(label_data_range, data.frame(
    band = b, x = x_pos, y = y_pos_range, label = label_range, stringsAsFactors = FALSE
  ))
}

# Plot sum_var fitting scatter plot
p_var_fit <- ggplot(grid_sf, aes(x = richness, y = sum_var)) +
  geom_point(aes(color = lat_band), alpha = 0.7, size = 2) +
  scale_color_manual(values = band_colors, name = "Latitude band") +
  scale_x_log10() +
  geom_smooth(aes(color = lat_band, group = lat_band), 
              method = "lm", se = TRUE, size = 0.8, alpha = 0.2) +
  geom_text(data = label_data_var,
            aes(x = x, y = y, label = label, color = band),
            size = 3.5, hjust = 1, vjust = 1, show.legend = FALSE) +
  labs(x = "Species richness per grid (log10 scale)", 
       y = "Sum of variances (sum_var)",
       title = "Morphological variance vs. species richness by latitudinal band") +
  theme_bw() + theme(legend.position = "right")

# Plot sum_range fitting scatter plot
p_range_fit <- ggplot(grid_sf, aes(x = richness, y = sum_range)) +
  geom_point(aes(color = lat_band), alpha = 0.7, size = 2) +
  scale_color_manual(values = band_colors, name = "Latitude band") +
  scale_x_log10() +
  geom_smooth(aes(color = lat_band, group = lat_band), 
              method = "lm", se = TRUE, size = 0.8, alpha = 0.2) +
  geom_text(data = label_data_range,
            aes(x = x, y = y, label = label, color = band),
            size = 3.5, hjust = 1, vjust = 1, show.legend = FALSE) +
  labs(x = "Species richness per grid (log10 scale)", 
       y = "Sum of ranges (sum_range)",
       title = "Morphological range vs. species richness by latitudinal band") +
  theme_bw() + theme(legend.position = "right")

print(p_var_fit)
print(p_range_fit)
ggsave("Richness_vs_SumVar_log10_linear_eq_icosa.pdf", p_var_fit, width = 10, height = 6)
ggsave("Richness_vs_SumRange_log10_linear_eq_icosa.pdf", p_range_fit, width = 10, height = 6)
