# Created by Xiaokang Liu on 24/09/2026 => xkliu@cug.edu.cn

packages <- c("data.table", "dplyr", "readxl", "ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
tax_path <- "../gastropods_taxonomy_combined.xlsx"
cat("读取分布数据...\n")
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
occ[, species := trimws(tolower(species))]


cat("读取分类表...\n")
taxonomy <- read_excel(tax_path)
taxonomy$species_clean <- trimws(tolower(taxonomy$species_clean))
taxonomy <- taxonomy[!duplicated(taxonomy$species_clean), ]

occ_tax <- merge(occ,
                 taxonomy[, c("species_clean", "Family")],
                 by.x = "species", by.y = "species_clean",
                 all.x = FALSE)
occ_tax <- occ_tax[!is.na(Family) & Family != ""]
cat("合并后记录数:", nrow(occ_tax), "\n")
cat("涉及科数:", uniqueN(occ_tax$Family), "\n")

occ_tax[, hemisphere := ifelse(latitude >= 0, "Northern", "Southern")]
occ_tax[, abs_lat    := abs(latitude)]

# 0-90°，每 5° 一个带
breaks <- seq(0, 90, by = 5)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])

occ_tax[, lat_band := cut(abs_lat,
                          breaks = breaks,
                          include.lowest = TRUE,
                          labels = labels)]

family_count_by_band <- occ_tax[, .(
  n_families = uniqueN(Family),
  n_species  = uniqueN(species),
  n_records  = .N
), by = .(hemisphere, lat_band)]

family_count_by_band <- family_count_by_band[order(hemisphere, lat_band)]

print(family_count_by_band)

p_family_lat <- ggplot(family_count_by_band,
                       aes(x = lat_band, y = n_families, fill = hemisphere)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_text(aes(label = n_families), vjust = -0.3, size = 3) +
  facet_wrap(~ hemisphere, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Northern" = "#e31a1c",
                               "Southern" = "#1f78b4")) +
  labs(x = "Absolute latitude band (degrees)",
       y = "Number of families",
       title = "Family richness along latitudinal gradient by hemisphere",
       subtitle = "Northern and Southern hemispheres plotted separately") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p_family_lat)
ggsave("family_richness_by_hemisphere.pdf", plot = p_family_lat,
       width = 10, height = 8, units = "in", device = "pdf")