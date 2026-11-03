# # Created by Xiaokang Liu on 17/08/2026 => xkliu@cug.edu.cn

library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(parallel)
library(viridis)


occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"

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

# ---------------------------- 2. grid cell --------------------------------
occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
bbox <- st_bbox(occ_sf)
cell_size <- 5

grid <- tryCatch({
  st_make_grid(occ_sf, cellsize = cell_size, square = FALSE, flat_topped = TRUE)
}, error = function(e) {
  st_make_grid(st_as_sfc(bbox), cellsize = cell_size, square = FALSE, flat_topped = TRUE)
})
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

sf_use_s2(FALSE)
valid <- st_is_valid(grid_sf)
if (!all(valid)) {
  grid_sf <- grid_sf[valid, ]
  grid_sf$grid_id <- 1:nrow(grid_sf)
}

intersections <- st_intersects(occ_sf, grid_sf)
occ$cell_id <- sapply(intersections, function(x) if (length(x) > 0) x[1] else NA)
occ <- occ[!is.na(cell_id)]

grid_summary <- occ[, .(
  occurrence_count = .N,
  species_list = list(unique(species))
), by = cell_id]
setnames(grid_summary, "cell_id", "grid_id")
grid_summary[, richness := lengths(species_list)]

grid_sf <- merge(grid_sf, grid_summary, by = "grid_id", all.x = TRUE)

# ---------------------------- 3. sum_var 和 sum_range --------------------
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

valid_idx <- which(lengths(grid_sf$species_list) >= 2)
grid_sf$sum_var <- NA
grid_sf$sum_range <- NA

if (length(valid_idx) > 100) {
  ncores <- max(1, detectCores() - 1)
  cat("Using", ncores, "cores for parallel computation...\n")
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("calc_disparity", "scores_10"))
    results <- parLapply(cl, grid_sf$species_list[valid_idx], function(sp) calc_disparity(sp, scores_10))
    stopCluster(cl)
  } else {
    results <- mclapply(grid_sf$species_list[valid_idx], function(sp) calc_disparity(sp, scores_10), mc.cores = ncores)
  }
  for (i in seq_along(valid_idx)) {
    grid_sf$sum_var[valid_idx[i]] <- results[[i]]["sum_var"]
    grid_sf$sum_range[valid_idx[i]] <- results[[i]]["sum_range"]
  }
} else {
  for (i in valid_idx) {
    res <- calc_disparity(grid_sf$species_list[[i]], scores_10)
    grid_sf$sum_var[i] <- res["sum_var"]
    grid_sf$sum_range[i] <- res["sum_range"]
  }
}
cat("Disparity calculation completed.\n")

# ---------------------------- 4. species > 10  ------------------------
grid_plot <- grid_sf[grid_sf$richness >= 10 & !is.na(grid_sf$richness), ]
cat("Grids with richness > 10:", nrow(grid_plot), "\n")


cat("\n=== Data ranges (positive values) ===\n")
cat("richness:  ", range(grid_plot$richness[grid_plot$richness > 0], na.rm = TRUE), "\n")
cat("sum_var:   ", range(grid_plot$sum_var[grid_plot$sum_var > 0], na.rm = TRUE), "\n")
cat("sum_range: ", range(grid_plot$sum_range[grid_plot$sum_range > 0], na.rm = TRUE), "\n")

# ---------------------------- 5. heatmaps -------------------------
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
      name = var_name
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = paste(title_prefix, var_name))
  return(p)
}


p_richness <- plot_heatmap(grid_plot, "richness")
p_sumvar   <- plot_heatmap(grid_plot, "sum_var")
p_sumrange <- plot_heatmap(grid_plot, "sum_range")


if (!is.null(p_richness)) {
  ggsave("Global_richness_heatmap.pdf", p_richness, width = 12, height = 6)
  print(p_richness)
}
if (!is.null(p_sumvar)) {
  ggsave("Global_sum_var_heatmap.pdf", p_sumvar, width = 12, height = 6)
  print(p_sumvar)
}
if (!is.null(p_sumrange)) {
  ggsave("Global_sum_range_heatmap.pdf", p_sumrange, width = 12, height = 6)
  print(p_sumrange)
}

cat("All heatmaps saved as PDF files.\n")
cat("Color mapping uses automatic data range (no limits).\n")



summary(grid_plot$richness)
summary(grid_plot$sum_var)
summary(grid_plot$sum_range)
