# Created by Xiaokang Liu on 11/10/2026 => xkliu@cug.edu.cn

library(geomorph)
library(dispRity)
library(readxl)
library(dplyr)

# ---------- Read PCA data and taxonomy ----------
PCA_data <- readRDS("./ModernGastropods/PCA_data.rds")
taxonomy <- read_excel("../gastropods_taxonomy_combined.xlsx")

scores <- PCA_data$x
stopifnot(nrow(scores) == nrow(taxonomy))

# ---------- Read IUCN data (marine sheet) ----------
iucn_path <- "../gastropoda all redlist_species_data summary/gastropods assessments in marine.xlsx"
iucn_raw <- read_excel(iucn_path, sheet = "marine")

col_lower <- tolower(colnames(iucn_raw))
sci_idx <- which(col_lower == "scientificname")
red_idx <- which(col_lower == "redlistcategory")
if (length(sci_idx) == 0) stop("scientificName column not found")
if (length(red_idx) == 0) stop("redlistCategory column not found")

iucn <- iucn_raw
colnames(iucn)[sci_idx] <- "scientificName"
colnames(iucn)[red_idx] <- "redlistCategory"

iucn <- iucn %>%
  mutate(scientificName = trimws(scientificName)) %>%
  filter(!is.na(scientificName), !is.na(redlistCategory))

group1 <- c("Least Concern")
group2 <- c("Near Threatened", "Lower Risk/near threatened", "Lower Risk/least concern")
group3 <- c("Critically Endangered", "Vulnerable", "Endangered", "Extinct")

iucn_grouped <- iucn %>%
  mutate(
    risk_group = case_when(
      redlistCategory %in% group1 ~ "Least Concern",
      redlistCategory %in% group2 ~ "Near Threatened / Lower Risk",
      redlistCategory %in% group3 ~ "Threatened / Extinct",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(risk_group)) %>%
  select(scientificName, risk_group)

cat("Number of valid IUCN grouped records:", nrow(iucn_grouped), "\n")

# ----------  Match risk groups to taxonomy ----------
species_col <- "species_clean"
if (!species_col %in% colnames(taxonomy)) {
  stop(paste("Column", species_col, "does not exist in taxonomy; available columns:",
             paste(colnames(taxonomy), collapse = ", ")))
}

taxonomy[[species_col]] <- trimws(taxonomy[[species_col]])
iucn_grouped$scientificName <- trimws(iucn_grouped$scientificName)

taxonomy$species_clean_lower <- tolower(taxonomy[[species_col]])
iucn_grouped$sci_lower <- tolower(iucn_grouped$scientificName)

taxonomy <- taxonomy %>%
  left_join(iucn_grouped, by = c("species_clean_lower" = "sci_lower"))

matched_count <- sum(!is.na(taxonomy$risk_group))
cat("Number of matched species (case-insensitive):", matched_count, "/", nrow(taxonomy), "\n")

taxonomy$risk_group[is.na(taxonomy$risk_group)] <- "Not Assessed"

risk_levels <- c("Least Concern", "Near Threatened / Lower Risk",
                 "Threatened / Extinct", "Not Assessed")
taxonomy$risk_group <- factor(taxonomy$risk_group, levels = risk_levels)

group_colors <- c("Least Concern" = "#2CA02C",
                  "Near Threatened / Lower Risk" = "#FFD700",
                  "Threatened / Extinct" = "#D62728",
                  "Not Assessed" = "gray70")

# ----------  Count species per family ----------
family_div <- taxonomy %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(!!sym(species_col))) %>%
  filter(n_species > 2) %>%
  pull(Family)

keep_family     <- taxonomy$Family %in% family_div
group_family    <- factor(taxonomy$Family[keep_family])
scores_family   <- scores[keep_family, ]
taxonomy_family <- taxonomy[keep_family, ]

# ---------- Compute normalized distances ----------
# For each family:
#   1) Use the mean of "all species" (including unassessed) on PC1-PC2 as the family centroid
#   2) Compute Euclidean distance from each species to the centroid
#   3) Normalize to 0-100% using the maximum distance within the family
compute_norm_dist <- function(scores_mat, family_vec, pc_dims = 1:10) {########################################
  n <- nrow(scores_mat)
  out <- rep(NA_real_, n)
  for (fam in unique(family_vec)) {
    idx <- which(family_vec == fam)
    if (length(idx) == 0) next
    if (length(idx) == 1) { out[idx] <- 0; next }
    
    pts <- scores_mat[idx, pc_dims, drop = FALSE]
    cen <- colMeans(pts)
    d   <- sqrt(rowSums((pts - matrix(cen, nrow(pts), ncol(pts), byrow = TRUE))^2))
    dmax <- max(d)
    out[idx] <- if (dmax > 0) d / dmax * 100 else 0
  }
  out
}

taxonomy_family$norm_dist_pct <- compute_norm_dist(scores_family, taxonomy_family$Family)

# ---------- 6. Single-family plotting function ----------
plot_family_norm_dist <- function(dist_pct, risk_vec, family_name, order_name,
                                  xlim = c(0, 100), ylim = c(0.5, 3.5)) {
  
  assessed_levels <- c("Least Concern",
                       "Near Threatened / Lower Risk",
                       "Threatened / Extinct")
  
  # Keep only assessed species
  keep <- !is.na(risk_vec) & as.character(risk_vec) %in% assessed_levels
  if (sum(keep) == 0) {
    cat("Skipping:", family_name, "(no assessed species)\n")
    return(invisible(NULL))
  }
  
  dist_use <- dist_pct[keep]
  risk_use <- droplevels(risk_vec[keep])
  
  pdf_file <- paste0("./ModernGastropods/IUCNfamily&ordermorphospace/",
                     "family_", order_name, "_", family_name,
                     "_normDist.pdf")
  pdf(pdf_file, width = 8, height = 7)
  
  # Vertical-axis "lanes"
  lane_map <- c("Least Concern"                 = 3,
                "Near Threatened / Lower Risk"  = 2,
                "Threatened / Extinct"          = 1)
  
  y_base <- lane_map[as.character(risk_use)]
  set.seed(42)
  y_jit  <- y_base + runif(length(y_base), -0.18, 0.18)
  
  # Empty plot
  plot(NA,
       xlim = xlim, ylim = ylim,
       xlab = "Normalized distance to family centroid (%)",
       ylab = "",
       yaxt = "n",
       main = paste0("FAMILY: ", family_name,
                     "\n(Order: ", order_name, ") - Assessed species only"))
  
  axis(2, at = 1:3,
       labels = c("Threatened / Extinct",
                  "Near Threatened / Lower Risk",
                  "Least Concern"),
       las = 1, cex.axis = 0.85)
  
  # Lane background lines
  abline(h = 1:3, col = "gray90", lty = 3)
  
  # Points
  for (lev in names(lane_map)) {
    idx <- which(as.character(risk_use) == lev)
    if (length(idx) > 0) {
      points(dist_use[idx], y_jit[idx],
             pch = 16, cex = 1.8,
             col = adjustcolor(group_colors[lev], alpha.f = 0.75))
    }
  }
  
  # Mean reference line for each lane
  for (lev in names(lane_map)) {
    idx <- which(as.character(risk_use) == lev)
    if (length(idx) > 0) {
      segments(mean(dist_use[idx]), lane_map[lev] - 0.28,
               mean(dist_use[idx]), lane_map[lev] + 0.28,
               col = group_colors[lev], lwd = 2, lty = 2)
    }
  }
  
  # Legend
  legend("topright",
         legend = c("Least Concern",
                    "Near Threatened / Lower Risk",
                    "Threatened / Extinct"),
         col = group_colors[c("Least Concern",
                              "Near Threatened / Lower Risk",
                              "Threatened / Extinct")],
         pch = 16, pt.cex = 1.5,
         title = "IUCN Status", bg = "white")
  
  n_table <- table(risk_use)
  legend("bottomright",
         legend = paste0(names(n_table), " = ", n_table),
         bty = "n", cex = 0.7)
  
  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# ---------- 7. Loop over each family ----------
dir.create("./ModernGastropods/IUCNfamily&ordermorphospace/",
           recursive = TRUE, showWarnings = FALSE)

for (fam in levels(group_family)) {
  order_vec <- unique(taxonomy_family$Order[taxonomy_family$Family == fam])
  if (length(order_vec) == 0) next
  order_name <- order_vec[1]
  
  idx <- which(taxonomy_family$Family == fam)
  dist_pct <- taxonomy_family$norm_dist_pct[idx]
  risk_vec <- taxonomy_family$risk_group[idx]
  
  plot_family_norm_dist(dist_pct, risk_vec, fam, order_name)
}

cat("Scatter plots of normalized distance for all families generated!\n")









# ============================================================
# [Summary plot] All assessed species from all families on one plot
#    X: normalized distance to family centroid (%)
#    Y: IUCN group (LC / NT-LR / Threatened-Extinct)
# ============================================================

assessed_levels <- c("Least Concern",
                     "Near Threatened / Lower Risk",
                     "Threatened / Extinct")

all_dist <- taxonomy_family$norm_dist_pct
all_risk <- as.character(taxonomy_family$risk_group)
all_family <- as.character(taxonomy_family$Family)

# Keep only assessed species
keep <- !is.na(all_dist) & all_risk %in% assessed_levels
dist_use <- all_dist[keep]
risk_use <- all_risk[keep]
fam_use  <- all_family[keep]

cat("Total points in summary plot:", length(dist_use), "\n")
print(table(risk_use))

# Lane mapping
lane_map <- c("Least Concern"                = 3,
              "Near Threatened / Lower Risk" = 2,
              "Threatened / Extinct"         = 1)

# Jitter (to reduce overlap along the same horizontal line)
set.seed(42)
y_base <- lane_map[risk_use]
y_jit  <- y_base + runif(length(y_base), -0.30, 0.30)

# ---------- Output PDF ----------
out_pdf <- "./ModernGastropods/IUCNfamily&ordermorphospace/ALL_families_normDist.pdf"
pdf(out_pdf, width = 10, height = 7)

par(mar = c(4.5, 8, 4, 2))

plot(NA,
     xlim = c(0, 100), ylim = c(0.5, 3.5),
     xlab = "Normalized distance to family centroid (%)",
     ylab = "", yaxt = "n",
     main = paste0("All families - Assessed species only  (n = ",
                   length(dist_use), ")"))

axis(2, at = 1:3,
     labels = c("Threatened / Extinct",
                "Near Threatened / Lower Risk",
                "Least Concern"),
     las = 1, cex.axis = 0.95)

abline(h = 1:3, col = "gray88", lty = 3)

# ---------- Points (semi-transparent, color = IUCN group) ----------
for (lev in names(lane_map)) {
  idx <- which(risk_use == lev)
  if (length(idx) > 0) {
    points(dist_use[idx], y_jit[idx],
           pch = 16, cex = 1.5,
           col = adjustcolor(group_colors[lev], alpha.f = 0.35))
  }
}

# ---------- Overlay statistics for each lane ----------
for (lev in names(lane_map)) {
  idx <- which(risk_use == lev)
  if (length(idx) == 0) next
  y0  <- lane_map[lev]
  col <- group_colors[lev]
  
  bp <- boxplot(dist_use[idx], plot = FALSE)
  
  # Median (thick solid line)
  segments(bp$stats[3], y0 - 0.34, bp$stats[3], y0 + 0.34,
           col = col, lwd = 3)
  # 25% / 75% quantiles (thin dashed lines)
  segments(bp$stats[2], y0 - 0.22, bp$stats[2], y0 + 0.22,
           col = col, lwd = 1.5, lty = 2)
  segments(bp$stats[4], y0 - 0.22, bp$stats[4], y0 + 0.22,
           col = col, lwd = 1.5, lty = 2)
  # Mean (open circle)
  points(mean(dist_use[idx]), y0,
         pch = 1, cex = 2, col = col, lwd = 2)
}

# ---------- Legend ----------
legend("topright",
       legend = c("Least Concern",
                  "Near Threatened / Lower Risk",
                  "Threatened / Extinct"),
       col = group_colors[c("Least Concern",
                            "Near Threatened / Lower Risk",
                            "Threatened / Extinct")],
       pch = 16, pt.cex = 1.5,
       title = "IUCN Status", bg = "white")

n_table <- table(risk_use)
legend("bottomright",
       legend = c(paste0(names(n_table), " = ", n_table),
                  "—  median",
                  "- -  25% / 75%",
                  "○  mean"),
       bty = "n", cex = 0.8)

dev.off()
cat("Saved:", out_pdf, "\n")




# ============================================================
# [Density plot] Probability density curves of normalized distance for the three IUCN groups
#    X: normalized distance to family centroid (%)
#    Y: relative probability density
# ============================================================

library(ggplot2)

assessed_levels <- c("Least Concern",
                     "Near Threatened / Lower Risk",
                     "Threatened / Extinct")

# Assemble plotting data (reusing the filtered results from section 8)
df_dens <- data.frame(
  norm_dist = dist_use,
  risk_group = factor(risk_use, levels = assessed_levels)
)

# Fixed color order (same as section 8)
col_map <- c("Least Concern"                 = "#2CA02C",
             "Near Threatened / Lower Risk"  = "#FFD700",
             "Threatened / Extinct"          = "#D62728")

# ---------- Density curves (with fill) ----------
p1 <- ggplot(df_dens, aes(x = norm_dist, fill = risk_group, color = risk_group)) +
  geom_density(alpha = 0.35, size = 1.1, adjust = 1) +
  scale_fill_manual(values  = col_map) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Probability density of normalized distance to family centroid",
    subtitle = paste0("Assessed species only (n = ", nrow(df_dens), ")"),
    x = "Normalized distance to family centroid (%)",
    y = "Relative probability density",
    fill = "IUCN Status", color = "IUCN Status"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave("./ModernGastropods/IUCNfamily&ordermorphospace/ALL_density_filled.pdf",
       p1, width = 8, height = 7)

# ---------- 9.2 Curves only (no fill, easier to compare peak positions) ----------
p2 <- ggplot(df_dens, aes(x = norm_dist, color = risk_group)) +
  geom_density(size = 1.2, adjust = 1) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Probability density of normalized distance to family centroid",
    subtitle = paste0("Assessed species only (n = ", nrow(df_dens), ")"),
    x = "Normalized distance to family centroid (%)",
    y = "Relative probability density",
    color = "IUCN Status"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.background = element_rect(fill = "white", color = "gray70"),
    plot.title = element_text(face = "bold")
  )

ggsave("./ModernGastropods/IUCNfamily&ordermorphospace/ALL_density_lines.pdf",
       p2, width = 8, height = 7)

cat("Density plots saved: ALL_density_filled.pdf / ALL_density_lines.pdf\n")
