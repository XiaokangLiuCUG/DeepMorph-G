# Created by Xiaokang Liu on 24/09/2026 => xkliu@cug.edu.cn

library(data.table)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(icosa)
library(sp)
library(iNEXT)

occ_path  <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"
out_dir    <- "./ModernGastropods/LDG"
if (!file.exists(occ_path))   stop("Occurrence file not found: ", occ_path)
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
  species_list     = list(unique(species))
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

# Keep grids with richness > 10 (consistent with previous code)
grid_sf <- grid_sf[!is.na(grid_sf$richness) & grid_sf$richness > 10, ]
cat("Number of retained grids (richness > 10):", nrow(grid_sf), "\n")

valid_grids <- grid_sf$grid_id


sp_records <- occ[cell_id %in% valid_grids,
                  .(n_rec = .N),
                  by = .(grid_id = cell_id, species)]

abund_list <- split(sp_records$n_rec, sp_records$grid_id)
names(abund_list) <- as.character(names(abund_list))
abund_list <- abund_list[sapply(abund_list, length) >= 2]
cat("Number of grids used for iNEXT computation:", length(abund_list), "\n")
cat("Running iNEXT for", length(abund_list), "grid cells ...\n")

inext_results <- lapply(seq_along(abund_list), function(i) {
  x <- abund_list[[i]]
  if (length(x) < 2) return(NULL)
  
  n_ref <- sum(x)    # reference sample size (total occurrences in this grid)
  
  res <- tryCatch({
    iNEXT(x,
          q        = 0,               # q = 0 -> species richness
          datatype = "abundance",     # abundance-based (individual-based)
          se       = FALSE,           # disable SE for speed
          nboot    = 0,
          endpoint = 3 * n_ref)       # * extrapolate to 3x reference sample size (Colwell's suggested upper bound)
  }, error = function(e) {
    message("Grid ", names(abund_list)[i], " failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(res)) return(NULL)
  
  dt <- as.data.table(res$iNextEst$size_based)
  dt[, grid_id := as.integer(names(abund_list)[i])]
  dt
})

rarefaction_dt <- rbindlist(inext_results, use.names = TRUE, fill = TRUE)
cat("Number of curve data rows:", nrow(rarefaction_dt), "\n")

# ---- Diagnostics: confirm the presence of Extrapolation ----
cat("Unique Method values: ",  paste(unique(rarefaction_dt$Method),  collapse = ", "), "\n")
cat("Unique Order.q values: ", paste(unique(rarefaction_dt$Order.q), collapse = ", "), "\n")

# ---- Filter: for q = 0, keep only Rarefaction / Extrapolation ----
rarefaction_dt <- rarefaction_dt[Order.q == 0]                     # * column name is Order.q
rarefaction_dt <- rarefaction_dt[Method %in% c("Rarefaction", "Extrapolation")]  
rarefaction_dt <- rarefaction_dt[!is.na(m) & !is.na(qD)]
cat("Number of rows after filtering:", nrow(rarefaction_dt), "\n")

grid_lat_dt <- grid_sf %>%
  mutate(lat_centroid = st_coordinates(st_point_on_surface(geometry))[, "Y"]) %>%
  st_set_geometry(NULL) %>%
  select(grid_id, lat_centroid) %>%
  as.data.table()

grid_lat_dt[, abs_centroid_lat := abs(lat_centroid)]

# Split into latitudinal bands using the specified breaks
grid_lat_dt[, lat_band := cut(abs_centroid_lat,
                              breaks        = c(0, 30, 45, 60, 90),
                              include.lowest = TRUE,
                              labels        = c("Lower", "Second", "Third", "Upper"))]

rarefaction_dt <- merge(rarefaction_dt,
                        grid_lat_dt[, .(grid_id, lat_band)],
                        by = "grid_id", all.x = TRUE)
rarefaction_dt <- rarefaction_dt[!is.na(lat_band)]

cat("Number of grids per latitudinal band:\n")
print(rarefaction_dt[, .(n_grids = uniqueN(grid_id)), by = lat_band])

# Export raw curve data
fwrite(rarefaction_dt,
       file.path(out_dir, "iNEXT_rarefaction_curves_raw.csv"))
cat("Saved:", file.path(out_dir, "iNEXT_rarefaction_curves_raw.csv"), "\n")

band_colors <- c("Lower"  = "#000080",
                 "Second" = "#33a02c",
                 "Third"  = "#ff7f00",
                 "Upper"  = "#6a3d9a")

rarefaction_dt[, lat_band := factor(lat_band,
                                    levels = c("Lower", "Second", "Third", "Upper"))]
rarefaction_dt[, Method   := factor(Method,
                                    levels = c("Rarefaction", "Extrapolation"))]

rarefaction_dt[, x_var := m]     # sample size
rarefaction_dt[, y_var := qD]    # q = 0 -> species richness


rarefaction_dt[, line_id := interaction(grid_id, Method, drop = TRUE)]

cat("Number of unique line_id values: ", uniqueN(rarefaction_dt$line_id), "\n")
cat("Expected ~ number of grids x 2 = ", uniqueN(rarefaction_dt$grid_id) * 2, "\n")

p_inext <- ggplot(rarefaction_dt,
                  aes(x        = x_var,
                      y        = y_var,
                      group    = line_id,      # * use line_id rather than grid_id
                      color    = lat_band,
                      linetype = Method)) +
  geom_line(alpha = 0.25, linewidth = 1) +
  # scale_x_log10(
  #   breaks = scales::log_breaks(n = 6),
  #   labels = scales::label_number(accuracy = 1)
  # ) +
  scale_color_manual(values = band_colors) +
  scale_linetype_manual(values = c("Rarefaction"  = "solid",
                                   "Extrapolation" = "22")) +
  facet_wrap(~ lat_band, ncol = 2, scales = "free") +
  labs(
    x = "Number of occurrences sampled (log10 scale)",
    y = "Expected species richness (q = 0)"
  ) +
  ggtitle("iNEXT rarefaction / extrapolation curves per grid cell across latitudinal bands",
          subtitle = "Abundance-based interpolation/extrapolation (q = 0); each thin line = one grid cell; solid = rarefaction, dashed = extrapolation") +
  theme_minimal(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line        = element_blank(),
    legend.position  = "none"
  )

ggsave(file.path(out_dir, "iNEXT_rarefaction_by_latband_waitui.pdf"),
       plot = p_inext, width = 11, height = 9, dpi = 300)

print(p_inext)

cat("\n============ finished ============\n")

cat("  ", file.path(out_dir, "iNEXT_rarefaction_curves_raw.csv"), "\n")
cat("  ", file.path(out_dir, "iNEXT_rarefaction_by_latband.pdf"),  "\n")