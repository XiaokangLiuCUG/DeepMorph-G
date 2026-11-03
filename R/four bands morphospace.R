# Created by Xiaokang Liu on 17/09/2026 => xkliu@cug.edu.cn

# ----------------------------- User-modified paths -----------------------------
morph_path <- "../ModernGastropods/PCA_data.rds"
occ_path <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
tax_path <- "../gastropods_taxonomy_combined.xlsx"

# Parameter settings
n_pcs <- 10

# --- 1. Extract the first two morphological PC column names (assumed PC1, PC2 or V1, V2) ---
# Since columns 5 onward of analysis_df are PCs, take the first two columns
pc1_col <- colnames(analysis_df)[5]
pc2_col <- colnames(analysis_df)[6]
cat("Morphological axes used:", pc1_col, "and", pc2_col, "\n")

# --- 2. Compute unified axis ranges (based on all species) ---
x_range <- range(analysis_df[[pc1_col]], na.rm = TRUE)
y_range <- range(analysis_df[[pc2_col]], na.rm = TRUE)
# Extend boundaries slightly (add 5% margin)
x_margin <- diff(x_range) * 0.05
y_margin <- diff(y_range) * 0.05
x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)

# --- 3. Compute convex hull for each latitudinal band ---
# Define function: compute convex hull polygon vertices
compute_hull <- function(df, x_col, y_col) {
  if (nrow(df) < 3) return(NULL)
  hull_idx <- chull(df[[x_col]], df[[y_col]])
  df[hull_idx, c(x_col, y_col)]
}

# Compute convex hulls grouped by latitudinal band
hull_data <- analysis_df %>%
  group_by(lat_band) %>%
  group_modify(~ compute_hull(.x, pc1_col, pc2_col)) %>%
  ungroup()

# If a band has too few species, it may be absent from hull_data; handle this
# Check which bands have no convex hull
present_bands <- unique(hull_data$lat_band)
missing_bands <- setdiff(levels(analysis_df$lat_band), present_bands)
if (length(missing_bands) > 0) {
  cat("Warning: the following bands have fewer than 3 species and cannot be plotted as convex hulls:", 
      paste(missing_bands, collapse = ", "), "\n")
}

# --- 4. Plot (facets, unified axes) ---
# Define colors for the four bands (same as before)
band_colors <- c("Lower"  = "#ff7f00",  
                 "Second" = "#33a02c", 
                 "Third"  = "#1f78b4",
                 "Upper"  = "#808080")

# Main plot: points + convex hull fill + convex hull border
p_morphospace <- ggplot(analysis_df, aes_string(x = pc1_col, y = pc2_col)) +
  # Plot points by latitudinal band (point transparency 0.3, colors distinguish bands)
  geom_point(aes(color = lat_band), alpha = 0.3, size = 1.5) +
  # Draw convex hulls (transparent fill, solid border)
  geom_polygon(data = hull_data, 
               aes_string(x = pc1_col, y = pc2_col, fill = "lat_band", color = "lat_band"),
               alpha = 0.2, size = 0.8, linetype = "solid") +
  # Set colors
  scale_color_manual(values = band_colors) +
  scale_fill_manual(values = band_colors) +
  # Facets: four bands arranged in rows
  facet_wrap(~ lat_band, ncol = 2) +
  # Unified axis ranges
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  # Labels and theme
  labs(x = "PC1", y = "PC2",
       title = "Morphospace occupation by latitudinal bands",
       subtitle = "Convex hulls show the occupied region; points show individual species") +
  theme_bw() +
  theme(legend.position = "none",  # Facets already label bands, so remove legend
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold", size = 12))

# Display plot
print(p_morphospace)


# --- 4. Plot (facets, unified axes) density plot ---
# Ensure ggpubr is loaded
library(ggpubr)

# Define colors
band_colors <- c("Lower"  = "#1f78b4", 
                 "Second" = "#33a02c",
                 "Third"  = "#ff7f00",
                 "Upper"  = "#e31a1c")

# Get latitudinal band order (assumed factor)
bands <- levels(analysis_df$lat_band)

# List to store subplots
plot_list <- list()

for (b in bands) {
  # Extract data for current band
  sub_df <- analysis_df[analysis_df$lat_band == b, ]
  
  # Compute convex hull (at least 3 points)
  if (nrow(sub_df) >= 3) {
    coords <- as.matrix(sub_df[, c(pc1_col, pc2_col)])
    hull_idx <- chull(coords)
    hull_coords <- coords[hull_idx, ]
    hull_df <- data.frame(x = hull_coords[, 1], y = hull_coords[, 2])
  } else {
    hull_df <- NULL
  }
  
  # Draw with ggscatterhist, specifying colors directly
  p <- ggscatterhist(
    data = sub_df,
    x = pc1_col, y = pc2_col,
    shape = 21,
    color = "black",                # Black point borders
    fill = band_colors[b],          # Fill directly with current band color
    size = 2,
    alpha = 0.6,
    margin.plot = "density",        # Marginal density plots
    margin.params = list(fill = band_colors[b], color = "black", size = 0.2),
    legend = "none",                # No legend for each subplot
    ggtheme = theme_bw(),
    xlab = "PC1", ylab = "PC2",
    title = b
  )
  
  # Add convex hull (if present)
  if (!is.null(hull_df)) {
    p <- p + 
      geom_polygon(data = hull_df, aes(x, y), 
                   fill = band_colors[b], alpha = 0.2, color = band_colors[b], 
                   linetype = "solid", size = 0.8, inherit.aes = FALSE)
  }
  
  plot_list[[b]] <- p
}

# Combine into a 2×2 grid with ggarrange (no shared legend)
p_combined <- ggarrange(plotlist = plot_list, 
                        ncol = 2, nrow = 2,
                        labels = "AUTO")

# Add global title and subtitle
p_combined <- annotate_figure(p_combined,
                              top = text_grob("Morphospace occupation by latitudinal bands\n(with marginal density)", 
                                              face = "bold", size = 16),
                              bottom = text_grob("Convex hulls show the occupied region; points show individual species", 
                                                 size = 10))

# Display plot
print(p_combined)




# --- 5. Optional: compute convex hull area for each band (morphospace size) ---
# Use the geometry package to compute area (requires installing geometry)
if (!require("geometry")) install.packages("geometry")
library(geometry)

hull_areas <- hull_data %>%
  group_by(lat_band) %>%
  summarise(area = convhulln(as.matrix(.[, c(pc1_col, pc2_col)]), "FA")$vol) %>%
  ungroup()
cat("\nConvex hull area for each band (morphospace size):\n")
print(hull_areas)

# --- 6. Save plot (optional) ---
ggsave("./ModernGastropods/Morphospace_by_lat_band.pdf", p_morphospace, width = 10, height = 8)




library(ggplot2)
library(dplyr)
# --- band 0-30, 30-45, 45-60, 60-90 ---
analysis_df <- analysis_df %>%
  mutate(
    lat_band = cut(
      mean_abs_lat,
      breaks = c(0, 30, 45, 60, 90),
      include.lowest = TRUE,
      labels = c("Lower", "Second", "Third", "Upper")
    )
  )

analysis_df$lat_band <- factor(
  analysis_df$lat_band,
  levels = c("Lower", "Second", "Third", "Upper")
)

# --- Unified colors: use the band colors you specified ---
band_colors <- c(
  "tropical"  = "#000080",   # dark blue
  "subtropical" = "#33a02c",   # green
  "temperate"  = "#ff7f00",   # orange
  "frigid"  = "#6a3d9a"    # purple
)

pc1_col <- colnames(analysis_df)[5]
pc2_col <- colnames(analysis_df)[6]
cat("Morphological axes used:", pc1_col, "and", pc2_col, "\n")

x_range <- range(analysis_df[[pc1_col]], na.rm = TRUE)
y_range <- range(analysis_df[[pc2_col]], na.rm = TRUE)

x_margin <- diff(x_range) * 0.05
y_margin <- diff(y_range) * 0.05

x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)

cat("Unified x-axis range:", x_limits, "\n")
cat("Unified y-axis range:", y_limits, "\n")


compute_hull <- function(df, x_col, y_col) {
  if (nrow(df) < 3) return(NULL)
  hull_idx <- chull(df[[x_col]], df[[y_col]])
  df[hull_idx, c(x_col, y_col)]
}

hull_data <- analysis_df %>%
  group_by(lat_band) %>%
  group_modify(~ compute_hull(.x, pc1_col, pc2_col)) %>%
  ungroup()

present_bands <- unique(hull_data$lat_band)
missing_bands <- setdiff(levels(analysis_df$lat_band), present_bands)

if (length(missing_bands) > 0) {
  cat("Warning: the following bands have fewer than 3 species and cannot be plotted as convex hulls:", 
      paste(missing_bands, collapse = ", "), "\n")
}

# ---  p_morphospace ---
p_morphospace <- ggplot(analysis_df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
  geom_point(aes(color = lat_band), alpha = 0.3, size = 1.5) +
  geom_polygon(
    data = hull_data,
    aes(x = .data[[pc1_col]], y = .data[[pc2_col]], fill = lat_band, color = lat_band),
    alpha = 0.2, size = 0.8, linetype = "solid"
  ) +
  scale_color_manual(values = band_colors) +
  scale_fill_manual(values = band_colors) +
  facet_wrap(~ lat_band, ncol = 2) +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  scale_x_continuous(breaks = pretty(x_limits, n = 5)) +
  scale_y_continuous(breaks = pretty(y_limits, n = 5)) +
  labs(
    x = "PC1", y = "PC2",
    title = "Morphospace occupation by latitudinal bands",
    subtitle = "Convex hulls show the occupied region; points show individual species"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_morphospace)

# --- save PDF ---
ggsave(
  filename = "./ModernGastropods/Morphospace_by_lat_band.pdf",
  plot = p_morphospace,
  width = 10,
  height = 8
)



# --- 6. Compute Sum of Range and Sum of Variance for the four latitudinal bands (first 10 PC axes) ---
library(dplyr)
library(tidyr)

n_pcs_use <- 10
pc_cols_10 <- colnames(analysis_df)[5:(4 + n_pcs_use)]

# Confirm columns exist
if (any(!pc_cols_10 %in% colnames(analysis_df))) {
  stop("The first 10 PC columns do not exist in analysis_df; please check the number of columns.")
}

band_stats <- analysis_df %>%
  group_by(lat_band) %>%
  summarise(
    across(
      all_of(pc_cols_10),
      list(
        range = ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE),
        var   = ~ var(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n_species = n(),
    .groups = "drop"
  )

sum_range <- band_stats %>%
  select(lat_band, n_species, ends_with("_range")) %>%
  mutate(sum_of_range = rowSums(across(ends_with("_range")), na.rm = TRUE)) %>%
  select(lat_band, n_species, sum_of_range)

sum_var <- band_stats %>%
  select(lat_band, ends_with("_var")) %>%
  mutate(sum_of_variance = rowSums(across(ends_with("_var")), na.rm = TRUE)) %>%
  select(lat_band, sum_of_variance)

summary_table <- sum_range %>%
  left_join(sum_var, by = "lat_band") %>%
  arrange(lat_band)

print(as.data.frame(summary_table), row.names = FALSE)

# save
write.csv(summary_table, "./ModernGastropods/lat_band_sum_range_variance.csv", row.names = FALSE)