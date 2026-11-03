# Created by Xiaokang Liu on 10/09/2026 => xkliu@cug.edu.cn
# ============================================================
# Modified version:
# 1) Density plot of species-to-family-centroid distances within each family
#    (normalized 0-1, using 0-100% for the density plot)
# 2) Proportion of species with distance > 0.8 vs species count within family
# ============================================================

library(geomorph)
library(dispRity)
library(readxl)
library(dplyr)
library(ggplot2)

# -------------------- 1. Data loading --------------------
PCA_data <- readRDS("./ModernGastropods/PCA_data.rds")
taxonomy <- read_excel("../gastropods_taxonomy_combined.xlsx")

scores <- PCA_data$x   # rows = species, columns = PCs
stopifnot(nrow(scores) == nrow(taxonomy))

# Number of PC axes used (first 10)
n_pc <- 10
if (ncol(scores) < n_pc) {
  warning("PCA scores have fewer than ", n_pc,
          " axes; using all ", ncol(scores), " axes.")
  n_pc <- ncol(scores)
}
pc_use <- 1:n_pc

# -------------------- 2. Filter families with > 10 species --------------------
family_stats <- taxonomy %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(species_clean), .groups = "drop") %>%
  filter(n_species > 10) %>%
  arrange(desc(n_species))

if (nrow(family_stats) == 0) {
  stop("No family with >10 species. Please lower the threshold.")
}

family_list <- family_stats$Family

# -------------------- 3. Compute species-to-family-centroid distances and normalize (first 10 PCs) --------------------
dist_list <- list()

for (fam in family_list) {
  idx <- which(taxonomy$Family == fam)
  
  # Take the first n_pc PCs (matrix: n_species x n_pc)
  pts <- scores[idx, pc_use, drop = FALSE]
  
  # Family centroid (n_pc-dimensional centroid)
  center <- colMeans(pts, na.rm = TRUE)
  
  # Multidimensional Euclidean distance: sqrt( sum_j (x_ij - c_j)^2 )
  d <- sqrt(rowSums(
    (pts - matrix(center, nrow(pts), ncol(pts), byrow = TRUE))^2,
    na.rm = TRUE))
  
  # Normalize to [0,1], then convert to percentage [0,100]
  rng <- range(d, na.rm = TRUE)
  if (diff(rng) > 0) {
    d_norm <- (d - rng[1]) / diff(rng)
  } else {
    d_norm <- rep(0, length(d))
  }
  d_pct <- d_norm * 100                    # <- new: percentage version
  
  dist_list[[fam]] <- data.frame(
    Family        = fam,
    species       = taxonomy$species_clean[idx],
    distance      = d,
    distance_norm = d_norm,
    distance_pct  = d_pct,                 # <- new column
    n_species     = length(idx),
    stringsAsFactors = FALSE
  )
}

dist_df <- do.call(rbind, dist_list)

# -------------------- 4. Density plot (normalized distance distribution, percentage x-axis) --------------------
# Colors: blue -> red -> yellow (low species count = blue, high species count = yellow)
p1 <- ggplot(dist_df,
             aes(x = distance_pct,                 # <- switched to percentage column
                 color = n_species,
                 group = Family)) +
  geom_density(size = 0.9, alpha = 0.25) +
  scale_color_gradientn(
    colours = c("blue", "red", "yellow"),
    name    = "Species\ncount"
  ) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +   # <- explicit tick marks
  labs(
    x = paste0("Normalized distance to family centroid (%, PC1\u2013PC", n_pc, ")"),
    y = "Density",
    title = paste0("Density of normalized distances (families with >10 species, PC1\u2013PC", n_pc, ")")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave("./ModernGastropods/Distance_density.pdf",
       p1, width = 9, height = 6)
cat("Saved: ./ModernGastropods/Distance_densitynew.pdf\n")

# -------------------- 5. Proportion with distance > 0.8 vs species count --------------------
ratio_df <- dist_df %>%
  group_by(Family, n_species) %>%
  summarise(
    n_far     = sum(distance_norm > 0.8),
    ratio_far = mean(distance_norm > 0.8),
    .groups   = "drop"
  )

# Binomial GLM: predicted values naturally fall in (0,1) and won't cross below y<0
fit <- glm(
  cbind(n_far, n_species - n_far) ~ n_species,
  data   = ratio_df,
  family = binomial(link = "logit")
)

# Smoothed prediction curve + 95% confidence band (response scale)
newdata <- data.frame(
  n_species = seq(min(ratio_df$n_species),
                  max(ratio_df$n_species),
                  length.out = 200)
)
pred_link <- predict(fit, newdata, type = "link", se.fit = TRUE)
inv_link  <- family(fit)$linkinv
newdata$fit <- inv_link(pred_link$fit)
newdata$lwr <- inv_link(pred_link$fit - 1.96 * pred_link$se.fit)
newdata$upr <- inv_link(pred_link$fit + 1.96 * pred_link$se.fit)

# Response-scale R²
r2 <- 1 - sum((ratio_df$ratio_far - fitted(fit))^2) /
  sum((ratio_df$ratio_far - mean(ratio_df$ratio_far))^2)

p2 <- ggplot(ratio_df, aes(x = n_species, y = ratio_far)) +
  geom_ribbon(data = newdata,
              aes(x = n_species, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              fill = "red", alpha = 0.15) +
  geom_line(data = newdata,
            aes(x = n_species, y = fit),
            inherit.aes = FALSE,
            color = "red", size = 1) +
  geom_point(size = 3, color = "steelblue") +
  annotate("text",
           x = max(ratio_df$n_species),
           y = max(ratio_df$ratio_far),
           label = paste0("R\u00B2 = ", round(r2, 3)),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(
    x = "Species count in family",
    y = "Proportion of species with normalized distance > 0.8",
    title = "Distance >0.8 proportion vs family species count",
    subtitle = paste0("Distances based on PC1\u2013PC", n_pc,
                      "; logistic regression")
  ) +
  theme_bw(base_size = 13)

ggsave("./ModernGastropods/Distance_rationew.pdf",
       p2, width = 8, height = 6)
cat("Saved: ./ModernGastropods/Distance_ratio.pdf\n")


# -------------------- 4b. Density plot (normalized distance distribution, log color scale) --------------------
# One density curve per family, colored by species count (log scale)
p1 <- ggplot(dist_df,
             aes(x = distance_pct,                 # <- switched to percentage column
                 color = n_species,
                 group = Family)) +
  geom_density(linewidth = 0.8, alpha = 0.25) +
  scale_color_gradientn(
    colours = c("blue", "red", "yellow"),
    trans   = "log10",
    name    = "Species\ncount"
  ) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +   # <- explicit tick marks
  labs(
    x = "Normalized distance to family centroid (%)",
    y = "Density",
    title = "Density of normalized distances (families with >10 species)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave("./ModernGastropods/Distance_densitylognew.pdf",
       p1, width = 9, height = 6)
cat("Saved: ./ModernGastropods/Distance_densitylog.pdf\n")


# -------------------- 5b. Proportion with distance > 0.8 vs species count (loess, x on log scale) --------------------
ratio_df <- dist_df %>%
  group_by(Family, n_species) %>%
  summarise(
    n_far     = sum(distance_norm > 0.8),
    ratio_far = mean(distance_norm > 0.8),
    .groups   = "drop"
  )

# Weighted loess: larger families have more reliable proportion estimates, so they get more weight
lo <- loess(ratio_far ~ n_species,
            data    = ratio_df,
            weights = n_species,
            span    = 0.75)

# Prediction grid changed to log-uniform to keep the curve smooth on the log axis
newdata <- data.frame(
  n_species = exp(seq(log(min(ratio_df$n_species)),
                      log(max(ratio_df$n_species)),
                      length.out = 300))
)
pred <- predict(lo, newdata, se = TRUE)

# Clamp to [0,1] to prevent the smoothed curve from going out of bounds
clamp01 <- function(x) pmin(pmax(x, 0), 1)
newdata$fit <- clamp01(pred$fit)
newdata$lwr <- clamp01(pred$fit - 1.96 * pred$se.fit)
newdata$upr <- clamp01(pred$fit + 1.96 * pred$se.fit)

# Response-scale R² (loess fitted values)
r2 <- 1 - sum((ratio_df$ratio_far - lo$fitted)^2) /
  sum((ratio_df$ratio_far - mean(ratio_df$ratio_far))^2)

x_r2 <- min(ratio_df$n_species)
y_r2 <- max(ratio_df$ratio_far)

p22 <- ggplot(ratio_df, aes(x = n_species, y = ratio_far)) +
  geom_ribbon(data = newdata,
              aes(x = n_species, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              fill = "red", alpha = 0.15) +
  geom_line(data = newdata,
            aes(x = n_species, y = fit),
            inherit.aes = FALSE,
            color = "red", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  scale_x_log10() +
  annotate("text", x = x_r2, y = y_r2,
           label = paste0("R\u00B2 = ", round(r2, 3)),
           hjust = -0.1, vjust = 1.2, size = 5) +
  labs(
    x = "Species count in family (log scale)",
    y = "Proportion of species with normalized distance > 0.8",
    title = "Distance >0.8 proportion vs family species count",
    subtitle = paste0("Distances based on PC1\u2013PC", n_pc,
                      "; weighted loess (span = 0.75)")
  ) +
  theme_bw(base_size = 13)

ggsave("./ModernGastropods/Distance_ratiologlooessnew.pdf",
       p22, width = 8, height = 6)
print(p22)


# -------------------- 5c. Proportion with distance > 0.8 vs species count (GLM, x on log scale) --------------------
ratio_df <- dist_df %>%
  group_by(Family, n_species) %>%
  summarise(
    n_far     = sum(distance_norm > 0.8),
    ratio_far = mean(distance_norm > 0.8),
    .groups   = "drop"
  )

# Binomial GLM
fit <- glm(
  cbind(n_far, n_species - n_far) ~ n_species,
  data   = ratio_df,
  family = binomial(link = "logit")
)

# Key point: prediction grid changed to log-uniform to keep the curve smooth on the log axis
newdata <- data.frame(
  n_species = exp(seq(log(min(ratio_df$n_species)),
                      log(max(ratio_df$n_species)),
                      length.out = 300))
)

pred_link <- predict(fit, newdata, type = "link", se.fit = TRUE)
inv_link  <- family(fit)$linkinv
newdata$fit <- inv_link(pred_link$fit)
newdata$lwr <- inv_link(pred_link$fit - 1.96 * pred_link$se.fit)
newdata$upr <- inv_link(pred_link$fit + 1.96 * pred_link$se.fit)

# Response-scale R²
r2 <- 1 - sum((ratio_df$ratio_far - fitted(fit))^2) /
  sum((ratio_df$ratio_far - mean(ratio_df$ratio_far))^2)

# On the log axis, choose a suitable position for R² (align to the upper-left using the minimum x)
x_r2 <- min(ratio_df$n_species)
y_r2 <- max(ratio_df$ratio_far)

p22 <- ggplot(ratio_df, aes(x = n_species, y = ratio_far)) +
  geom_ribbon(data = newdata,
              aes(x = n_species, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              fill = "red", alpha = 0.15) +
  geom_line(data = newdata,
            aes(x = n_species, y = fit),
            inherit.aes = FALSE,
            color = "red", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  scale_x_log10() +                       # <- key point: x-axis on log scale
  annotate("text",
           x = x_r2, y = y_r2,
           label = paste0("R\u00B2 = ", round(r2, 3)),
           hjust = -0.1, vjust = 1.2, size = 5) +
  labs(
    x = "Species count in family (log scale)",
    y = "Proportion of species with normalized distance > 0.8",
    title = "Distance >0.8 proportion vs family species count",
    subtitle = paste0("Distances based on PC1\u2013PC", n_pc,
                      "; logistic regression (log-x axis)")
  ) +
  theme_bw(base_size = 13)

ggsave("./ModernGastropods/Distance_ratiolognew.pdf",
       p22, width = 8, height = 6)
print(p22)
cat("Saved: ./ModernGastropods/Distance_ratiolog.pdf\n")