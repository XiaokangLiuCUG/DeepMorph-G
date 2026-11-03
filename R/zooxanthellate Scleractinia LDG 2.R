# Created by Xiaokang Liu on 30/06/2026 => xkliu@cug.edu.cn
# Species richness of reef-building corals along a latitudinal gradient
# Load required packages
packages <- c("data.table", "ggplot2", "dplyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ----------------------------- File paths -----------------------------
occ_file <- "../Scleractinia 0020924-260721160103020/0020924-260721160103020.csv"
species_list_file <- "../zooxanthellate Scleractinia species.txt"

# ----------------------------- 1. Read the reef-building coral species list -----------------------------
cat("Reading reef-building coral species list...\n")
species_raw <- readLines(species_list_file, warn = FALSE)
species_raw <- trimws(species_raw)
species_raw <- species_raw[species_raw != ""]
reef_species <- data.table(species_raw = species_raw)
reef_species[, species := tolower(species_raw)]
reef_species <- unique(reef_species, by = "species")
cat("Number of reef-building coral species:", nrow(reef_species), "\n")

# ----------------------------- 2. Read occurrence records (including depth) -----------------------------
cat("Reading occurrence file (large file, ~8 GB, please be patient)...\n")
occ <- fread(
  occ_file,
  select = c("species", "decimalLatitude", "decimalLongitude", "depth"),
  na.strings = c("", "NA", "NULL"),
  stringsAsFactors = FALSE
)
cat("Number of raw records:", nrow(occ), "\n")

# Coordinate cleaning: remove records with missing coordinates or species names
occ <- occ[complete.cases(occ[, .(species, decimalLatitude, decimalLongitude)])]
occ <- occ[decimalLatitude >= -90 & decimalLatitude <= 90 & 
             decimalLongitude >= -180 & decimalLongitude <= 180]
occ[, species := trimws(tolower(species))]
cat("Number of records after coordinate cleaning:", nrow(occ), "\n")

# ----------------------------- 3. Depth filtering (keep depth <= 200 or missing) -----------------------------
cat("Filtering records with depth > 200...\n")
occ <- occ[is.na(depth) | depth <= 200]
cat("Number of records after depth filtering:", nrow(occ), "\n")

# ----------------------------- 4. Match reef-building coral species -----------------------------
cat("Matching reef-building coral species...\n")
setkey(occ, species)
setkey(reef_species, species)
reef_occ <- occ[reef_species, nomatch = 0]  # keep only matched records
cat("Number of records after matching:", nrow(reef_occ), "\n")

rm(occ); gc()

# Compute mean and standard deviation of absolute latitude for each species
species_stats <- reef_occ[, .(mean_abs_lat = mean(abs(decimalLatitude), na.rm = TRUE),
                              sd_abs_lat = sd(abs(decimalLatitude), na.rm = TRUE)),
                          by = species]

# Merge back into the original data
reef_occ <- merge(reef_occ, species_stats, by = "species")

# Filter condition: keep if |lat| ≤ 30, or if |lat| is within mean_abs ± 3*sd_abs
reef_occ <- reef_occ[abs(decimalLatitude) <= 30 | 
                       (abs(decimalLatitude) >= (mean_abs_lat - 3 * sd_abs_lat) &
                          abs(decimalLatitude) <= (mean_abs_lat + 3 * sd_abs_lat))]
cat("Number of records after outlier removal:", nrow(reef_occ), "\n")

# Optional: remove species with sd_abs_lat = NA (only 1-2 records); uncomment to apply
# reef_occ <- reef_occ[!is.na(sd_abs_lat)]

# Remove helper columns
reef_occ[, c("mean_abs_lat", "sd_abs_lat") := NULL]
# ==============================================================================

# ----------------------------- 5. Define latitude bins (4-degree resolution) -----------------------------
resolution <- 4
reef_occ[, lat_bin := floor(decimalLatitude / resolution) * resolution]

# ----------------------------- 6. Compute species richness per latitude bin (unique species count) -----------------------------
richness <- reef_occ[, .(n_species = uniqueN(species)), by = lat_bin]
setorder(richness, lat_bin)
cat("Number of latitude bins:", nrow(richness), "\n")

# ----------------------------- 7. Plot -----------------------------
p <- ggplot(richness, aes(x = lat_bin, y = n_species)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = resolution * 0.8) +
  scale_x_continuous(breaks = seq(-90, 90, by = 10), limits = c(-90, 90)) +
  labs(
    title = "Species richness of zooxanthellate Scleractinia (|lat|≤30° protected)",
    subtitle = paste("Depth ≤ 200 m | Resolution:", resolution, "° | ",
                     "Removed records with |lat| > 30° beyond mean_abs ± 3*sd_abs per species"),
    x = "Latitude (°)",
    y = "Number of species"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p)

# Optional save
ggsave("reef_coral_richness_abs_lat_outlier_removed_protect_lowlat.pdf", p, width = 10, height = 6)

cat("\nDone! Figure displayed and saved as PDF.\n")




# ----------------------------- 5. Define latitude bins (4-degree resolution) -----------------------------
resolution <- 4
reef_occ[, lat_bin := floor(decimalLatitude / resolution) * resolution]

# ----------------------------- 6. Compute species richness per latitude bin (unique species count) -----------------------------
richness <- reef_occ[, .(n_species = uniqueN(species)), by = lat_bin]
setorder(richness, lat_bin)
cat("Number of latitude bins:", nrow(richness), "\n")

# ----------------------------- 7. Plot (bar chart + smoothed trend line) -----------------------------
p <- ggplot(richness, aes(x = lat_bin, y = n_species)) +
  geom_col(fill = "steelblue", alpha = 0.6, width = resolution) +# * 0.8
  geom_smooth(method = "loess", span = 0.5, se = TRUE, color = "red", size = 1.2) +
  scale_x_continuous(breaks = seq(-90, 90, by = 10), limits = c(-90, 90)) +
  scale_y_continuous(limits = c(0, NA)) +  # force y-axis to start at 0
  labs(
    title = "Species richness of zooxanthellate Scleractinia (|lat|≤30° protected)",
    subtitle = paste("Depth ≤ 200 m | Resolution:", resolution, "° | LOESS (span=0.7)"),
    x = "Latitude (°)",
    y = "Number of species"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p)

# Optional save
ggsave("reef_coral_richness_abs_lat_outlier_removed_protect_lowlat.pdf", p, width = 10, height = 6)############################Use
