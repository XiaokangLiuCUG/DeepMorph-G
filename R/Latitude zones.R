# Created by Xiaokang Liu on 02/08/2026 => xkliu@cug.edu.cn

library(data.table)
occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
cat("Total valid occurrence records:", nrow(occ), "\n")
if (nrow(occ) == 0) stop("No valid coordinates after filtering.")


species_stats <- occ[, .(
  min_lat      = min(latitude),
  max_lat      = max(latitude),
  min_abs_lat  = min(abs(latitude)),
  max_abs_lat  = max(abs(latitude)),
  mean_abs_lat = mean(abs(latitude)),
  n_records    = .N
), by = species]


species_stats[, lat_band := cut(mean_abs_lat,
                                breaks = c(0, 30, 45, 60, 90),
                                include.lowest = TRUE,
                                labels = c("Lower", "Second", "Third", "Upper"))]
setorder(species_stats, mean_abs_lat)
morph_path <- "D:/OneDrive/文档/ModernGastropods/PCA_data.rds"
if (!file.exists(morph_path)) stop("PCA_data.rds file does not exist; please check the path.")
pca_obj <- readRDS(morph_path)

if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else if (is.matrix(pca_obj) || is.data.frame(pca_obj)) {
  morph <- as.matrix(pca_obj)
} else {
  stop("PCA_data.rds format not recognized; please confirm it is a prcomp or matrix object.")
}
if (is.null(rownames(morph))) stop("Morphological matrix rownames are empty; cannot match species.")
cat("Morphological matrix dimensions:", dim(morph), "\n")

# ----------------------------- 5. Match latitudinal bands -----------------------------
band_df <- species_stats[, .(species, lat_band)]
all_morph_species <- rownames(morph)

morph_bands <- band_df[match(all_morph_species, band_df$species), lat_band]
morph_bands <- as.character(morph_bands)
morph_bands[is.na(morph_bands)] <- "unassigned"

# Check matching status
cat("\nNumber of species per band (including unassigned):\n")
print(table(morph_bands))

# Extract the first two principal components
scores <- morph[, 1:2, drop = FALSE]

band_colors <- c("Lower"      = "#000080", 
                 "Second"     = "#33a02c", 
                 "Third"      = "#ff7f00", 
                 "Upper"      = "#6a3d9a",
                 "unassigned" = "#666666")
pdf("Morphospace_by_LatBand_ordered.pdf", width = 10, height = 7.5)##########################################################
x_range <- range(scores[,1], na.rm = TRUE)
y_range <- range(scores[,2], na.rm = TRUE)
plot(NA, xlim = x_range, ylim = y_range,
     xlab = "PC1", ylab = "PC2",
     main = "Morphospace (PC1 vs PC2) colored by latitudinal band")


unassigned_idx <- morph_bands == "unassigned"
if (any(unassigned_idx)) {
  points(scores[unassigned_idx, 1], scores[unassigned_idx, 2],
         col = band_colors["unassigned"], pch = 16, cex = 0.7)
}

order_bands <- c("Lower", "Second", "Third", "Upper")
for (b in order_bands) {
  idx <- morph_bands == b
  if (any(idx)) {
    points(scores[idx, 1], scores[idx, 2],
           col = band_colors[b], pch = 16, cex = 0.7)
  }
}

legend("topright", 
       legend = c("Lower", "Second", "Third", "Upper", "unassigned"), 
       col = band_colors[c("Lower", "Second", "Third", "Upper", "unassigned")], 
       pch = 19, cex = 1.2, 
       title = "Latitude band")

dev.off()
cat("Morphospace plot saved to Morphospace_by_LatBand_ordered.pdf\n")


# New code
band_counts <- table(factor(morph_bands, levels = names(band_colors)))
cat("\nNumber of species per band (ordered by band_colors):\n")
print(band_counts)

print(data.frame(
  band = names(band_counts),
  n = as.integer(band_counts)
))


# Load required packages
packages <- c("data.table", "dplyr", "readxl", "ggplot2", "ggrepel", "nlme")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


morph_path <- "../ModernGastropods/PCA_data.rds"
occ_path <- "..Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
tax_path <- "../gastropods_taxonomy_combined.xlsx"

# Parameter settings
n_pcs <- 10
min_species_per_family <- 10

cat("Reading distribution data...\n")
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude", "Family", "Order"))
setnames(occ, c("species", "latitude", "longitude", "family", "order"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
occ[, species := trimws(tolower(species))]
occ[, family := trimws(family)]
occ[, order := trimws(order)]

# Summarize by species
species_stats <- occ[, .(
  mean_abs_lat = mean(abs(latitude)),
  n_records = .N,
  family = first(na.omit(family)),
  order = first(na.omit(order))
), by = species]


species_stats[, lat_band := cut(mean_abs_lat,
                                breaks = c(0, 30, 45, 60, 90),
                                include.lowest = TRUE,
                                labels = c("Lower", "Second", "Third", "Upper"))]


cat("Reading morphological data...\n")
pca_obj <- readRDS(morph_path)
if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else {
  morph <- as.matrix(pca_obj)
}
rownames(morph) <- trimws(tolower(rownames(morph)))
cat("Number of species in morphological matrix:", nrow(morph), "\n")

cat("Reading taxonomy table...\n")
taxonomy <- read_excel(tax_path)
taxonomy$species_clean <- trimws(tolower(taxonomy$species_clean))
taxonomy <- taxonomy[!duplicated(taxonomy$species_clean), ]

# ----------------------------- Take the intersection of all three -----------------------------
common_species <- intersect(rownames(morph), species_stats$species)
common_species <- intersect(common_species, taxonomy$species_clean)
cat("Number of species shared by all three datasets:", length(common_species), "\n")
if (length(common_species) == 0) stop("No shared species.")

morph_common <- morph[common_species, ]
species_stats_common <- species_stats[match(common_species, species_stats$species), ]

analysis_df <- data.frame(
  species = common_species,
  family = species_stats_common$family,
  order = species_stats_common$order,
  mean_abs_lat = species_stats_common$mean_abs_lat,
  lat_band = species_stats_common$lat_band,
  morph_common[, 1:n_pcs, drop = FALSE],
  stringsAsFactors = FALSE
)
analysis_df <- analysis_df[!is.na(analysis_df$family) & !is.na(analysis_df$order), ]

# ----------------------------- Filter by family -----------------------------
family_counts <- analysis_df %>%
  group_by(family) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > min_species_per_family)

families_ok <- family_counts$family
analysis_df <- analysis_df[analysis_df$family %in% families_ok, ]
cat("Number of species after filtering:", nrow(analysis_df), "\n")
cat("Number of qualifying families:", length(families_ok), "\n")
if (nrow(analysis_df) == 0) stop("No qualifying families.")


pc_cols <- 6:(5 + n_pcs)
global_centroid <- colMeans(analysis_df[, pc_cols, drop = FALSE], na.rm = TRUE)
pc_mat_all <- as.matrix(analysis_df[, pc_cols])
analysis_df$dist_to_centroid <- sqrt(rowSums(sweep(pc_mat_all, 2, global_centroid)^2))

cor_pearson <- cor.test(analysis_df$dist_to_centroid, analysis_df$mean_abs_lat, method = "pearson")
cor_spearman <- cor.test(analysis_df$dist_to_centroid, analysis_df$mean_abs_lat, method = "spearman", exact = FALSE)

cat("\n=== Pearson correlation (global centroid) ===\n")
print(cor_pearson)
cat("\n=== Spearman rank correlation (global centroid) ===\n")
print(cor_spearman)


band_colors <- c("Lower"   = "#000080",   # dark blue
                 "Second"  = "#33a02c",   # green
                 "Third"   = "#ff7f00",   # orange
                 "Upper"   = "#6a3d9a")   # purple

# ----------------------------- Scatter plot -----------------------------
p_band <- ggplot(analysis_df, aes(x = mean_abs_lat, y = dist_to_centroid, color = lat_band)) +
  geom_point(alpha = 0.6, size = 1, stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed",
              aes(group = 1), formula = y ~ x, size = 1) +
  scale_color_manual(values = band_colors, name = "Latitude band") +
  labs(x = "Mean absolute latitude (degrees)",
       y = "Euclidean distance to global centroid",
       title = "Morphological distinctiveness vs. mean absolute latitude",
       subtitle = paste0("Pearson r = ", round(cor_pearson$estimate, 3),
                         ", p = ", format(cor_pearson$p.value, digits = 3))) +
  theme_bw() +
  theme(legend.position = "right")
print(p_band)

# ----------------------------- Boxplot -----------------------------
p_box <- ggplot(analysis_df, aes(x = lat_band, y = dist_to_centroid, fill = lat_band)) +
  geom_boxplot() +
  scale_fill_manual(values = band_colors) +
  labs(x = "Latitude band", y = "Distance to centroid",
       title = "Morphological distance distribution by latitude band") +
  theme_bw() +
  theme(legend.position = "none")
print(p_box)


ggsave("Scatter_Distance_vs_LatBand.pdf", p_band, width = 10, height = 6)##########################################################
ggsave("Boxplot_Distance_by_LatBand.pdf", p_box, width = 10, height = 6)##########################################################







library(data.table)
# 1. Read occurrence data ------------------------------------------------
occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))

# Filter valid records
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
cat("Total valid occurrence records:", nrow(occ), "\n")
if (nrow(occ) == 0) stop("No valid coordinates after filtering.")

# 2. Summarize by species: compute mean_abs_lat for each species -----------------------------
species_stats <- occ[, .(
  min_lat      = min(latitude),
  max_lat      = max(latitude),
  min_abs_lat  = min(abs(latitude)),
  max_abs_lat  = max(abs(latitude)),
  mean_abs_lat = mean(abs(latitude)),   # core metric: mean absolute latitude of all records
  n_records    = .N
), by = species]

species_stats[, lat_band := cut(mean_abs_lat,
                                breaks = c(0, 30, 45, 60, 90),
                                include.lowest = TRUE,
                                labels = c("Lower", "Second", "Third", "Upper")
)]

setorder(species_stats, mean_abs_lat)
band_summary <- species_stats[, .(
  n_species = .N,
  min_mean_abs = min(mean_abs_lat),
  max_mean_abs = max(mean_abs_lat)
), by = lat_band]
print(band_summary)
fwrite(species_stats, "species_lat_bands_MeanAbs.csv", bom = FALSE)




# ========== Visualization of latitudinal bands in morphospace ==========
morph_path <- "D:/OneDrive/文档/ModernGastropods/PCA_data.rds"
if (!file.exists(morph_path)) stop("PCA_data.rds file does not exist; please check the path.")
pca_obj <- readRDS(morph_path)

if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else if (is.matrix(pca_obj) || is.data.frame(pca_obj)) {
  morph <- as.matrix(pca_obj)
} else {
  stop("PCA_data.rds format not recognized; please confirm it is a prcomp or matrix object.")
}
if (is.null(rownames(morph))) stop("Morphological matrix rownames are empty; cannot match species.")
cat("Morphological matrix dimensions:", dim(morph), "\n")
band_df <- species_stats[, .(species, lat_band)]

common_species <- intersect(rownames(morph), band_df$species)
if (length(common_species) == 0) stop("No common species between morphology and occurrence.")
cat("Number of species with morphology data:", length(common_species), "\n")

# Extract scores on the first two principal components (PC1, PC2)
scores <- morph[common_species, 1:2, drop = FALSE]
# Get the corresponding latitudinal band labels in matching order
bands <- band_df[match(common_species, band_df$species), lat_band]


# ========== Latitudinal distribution range vs morphological distance ==========

morph_common <- morph[common_species, , drop = FALSE]  # rownames = species, columns = PC1, PC2, ...

centroid <- colMeans(morph_common)

dist_to_centroid <- sqrt(rowSums((morph_common - matrix(centroid,
                                                        nrow = nrow(morph_common),
                                                        ncol = ncol(morph_common),
                                                        byrow = TRUE))^2))

plot_df <- species_stats[match(common_species, species_stats$species),
                         .(species, min_lat, max_lat, lat_band)]
plot_df$distance <- dist_to_centroid
# Sort by distance in descending order (larger distances on top)
plot_df <- plot_df[order(-distance), ]


band_colors <- c("Lower"   = "#000080",   # dark blue
                 "Second"  = "#33a02c",   # green
                 "Third"   = "#ff7f00",   # orange
                 "Upper"   = "#6a3d9a")   # purple
# 6. Plot and save as PDF
pdf("Latitudinal_Range_vs_Morph_Distance.pdf", width = 10, height = 10)##########################################################


plot(NA, xlim = c(-90, 90),
     ylim = c(0, max(plot_df$distance) * 1.05),
     xlab = "Latitude (°)",
     ylab = "Morphological distance to centroid",
     main = "Latitudinal range vs morphological distance",
     xaxt = "n")  # hide default x-axis

# Custom x-axis ticks (-90, -60, -30, 0, 30, 60, 90)
axis(side = 1, at = seq(-90, 90, by = 30))

# Draw line segments (existing loop unchanged)
for (i in 1:nrow(plot_df)) {
  segments(x0 = plot_df$min_lat[i],
           x1 = plot_df$max_lat[i],
           y0 = plot_df$distance[i],
           y1 = plot_df$distance[i],
           col = band_colors[plot_df$lat_band[i]],
           lwd = 0.8)
}
legend("topright", legend = names(band_colors), col = band_colors, lwd = 2, title = "Latitude band")
dev.off()