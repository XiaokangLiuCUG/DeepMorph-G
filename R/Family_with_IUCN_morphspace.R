# Created by Xiaokang Liu on 23/07/2026 => xkliu@cug.edu.cn

library(geomorph)
library(dispRity)
library(readxl)
library(dplyr)

PCA_data <- readRDS("./ModernGastropods/PCA_data.rds")
taxonomy <- read_excel("../gastropods_taxonomy_combined.xlsx")

scores <- PCA_data$x   # rows = species, columns = PCs
stopifnot(nrow(scores) == nrow(taxonomy))

iucn_path <- "../gastropoda all redlist_species_data summary/assessments marine.xlsx"
iucn_raw <- read_excel(iucn_path, sheet = "marine")

# Automatically adapt column names (case-insensitive)
col_lower <- tolower(colnames(iucn_raw))
sci_idx <- which(col_lower == "scientificname")
red_idx <- which(col_lower == "redlistcategory")
if (length(sci_idx) == 0) stop("scientificName column not found")
if (length(red_idx) == 0) stop("redlistCategory column not found")

iucn <- iucn_raw
colnames(iucn)[sci_idx] <- "scientificName"
colnames(iucn)[red_idx] <- "redlistCategory"

# Trim whitespace and remove missing values
iucn <- iucn %>%
  mutate(scientificName = trimws(scientificName)) %>%
  filter(!is.na(scientificName), !is.na(redlistCategory))

# Define the three grouping categories
group1 <- c("Least Concern")
group2 <- c("Near Threatened", "Lower Risk/near threatened", "Lower Risk/least concern")
group3 <- c("Critically Endangered", "Vulnerable", "Endangered", "Extinct")

# Create risk group mapping (data frame: species, risk_group)
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

species_col <- "species_clean"  
if (!species_col %in% colnames(taxonomy)) {
  stop(paste("Column", species_col, "does not exist in taxonomy; available columns: ", 
             paste(colnames(taxonomy), collapse = ", ")))
}

taxonomy[[species_col]] <- trimws(taxonomy[[species_col]])
iucn_grouped$scientificName <- trimws(iucn_grouped$scientificName)
taxonomy$species_clean_lower <- tolower(taxonomy[[species_col]])
iucn_grouped$sci_lower <- tolower(iucn_grouped$scientificName)
taxonomy <- taxonomy %>%
  left_join(iucn_grouped, by = c("species_clean_lower" = "sci_lower"))

matched_count <- sum(!is.na(taxonomy$risk_group))
cat("\n=== Number of matched species (case-insensitive) ===\n", matched_count, "/", nrow(taxonomy), "\n")

# Set unmatched entries to "Not Assessed"
taxonomy$risk_group[is.na(taxonomy$risk_group)] <- "Not Assessed"

# Convert risk_group to a factor with a fixed order
risk_levels <- c("Least Concern", "Near Threatened / Lower Risk", 
                 "Threatened / Extinct", "Not Assessed")
taxonomy$risk_group <- factor(taxonomy$risk_group, levels = risk_levels)

# Define color mapping (three categories + gray)
group_colors <- c("Least Concern" = "#2CA02C",
                  "Near Threatened / Lower Risk" = "#FFD700",
                  "Threatened / Extinct" = "#D62728",
                  "Not Assessed" = "gray70")

# ---------- Count the number of species in each family ----------
family_div <- taxonomy %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(!!sym(species_col))) %>%
  filter(n_species > 2) %>%
  pull(Family)

# Filter retained families
keep_family <- taxonomy$Family %in% family_div
group_family <- factor(taxonomy$Family[keep_family])
scores_family <- scores[keep_family, ]
taxonomy_family <- taxonomy[keep_family, ]

# ---------- Define plotting function (single family) ----------
plot_family_group <- function(pc1, pc2, risk_vec, family_name, order_name,
                              xlim = c(-0.6, 0.55), ylim = c(-0.7, 0.5)) {
  
  pdf_file <- paste0("./ModernGastropods/IUCNfamily&ordermorphospace/family_", 
                     order_name, "_", family_name, ".pdf")
  pdf(pdf_file, width = 8, height = 6)
  plot(0, 0, type = "n", 
       xlab = "PC1", ylab = "PC2",
       xlim = xlim, ylim = ylim,
       main = paste0("FAMILY: ", family_name, "\n(Order: ", order_name, ")"))

  for (lev in levels(risk_vec)) {
    idx <- which(risk_vec == lev)
    if (length(idx) > 0) {
      points(pc1[idx], pc2[idx], 
             pch = 0, cex = 1.5, 
             col = group_colors[lev])
    }
  }
  
  # Plot the centroid (based on all points, red filled square)
  points(mean(pc1), mean(pc2), pch = 15, cex = 2, col = "red")
  
  # Plot the convex hull (based on all points, requires at least 3 points)
  if (length(pc1) >= 3) {
    hpts <- chull(pc1, pc2)
    hpts <- c(hpts, hpts[1])
    lines(pc1[hpts], pc2[hpts], col = "blue", lwd = 2)
    polygon(pc1[hpts], pc2[hpts], col = adjustcolor("blue", alpha.f = 0.2), border = NA)
  } else {
    text(mean(pc1), mean(pc2), "n < 3\n(no hull)", cex = 0.8, col = "gray40")
  }
  
  legend("topright", 
         legend = c("Least Concern", "Near Threatened / Lower Risk", "Threatened / Extinct"),
         col = group_colors[c("Least Concern", "Near Threatened / Lower Risk", "Threatened / Extinct")],
         pch = 0, pt.cex = 1.5, title = "IUCN status")
  
  n_table <- table(risk_vec)
  n_text <- paste0(names(n_table), " = ", n_table)
  legend("bottomright", legend = n_text, bty = "n", cex = 0.7)
  
  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# ---------- Loop over each family ----------
dir.create("./ModernGastropods/IUCNfamily&ordermorphospace/", recursive = TRUE, showWarnings = FALSE)

family_names <- levels(group_family)
for (fam in family_names) {
  # Get the Order to which this Family belongs
  order_vec <- unique(taxonomy_family$Order[taxonomy_family$Family == fam])
  if (length(order_vec) == 0) next
  order_name <- order_vec[1]
  
  # Extract the data indices for this family (in taxonomy_family)
  idx <- which(taxonomy_family$Family == fam)
  pc1 <- scores_family[idx, 1]
  pc2 <- scores_family[idx, 2]
  risk_vec <- taxonomy_family$risk_group[idx]   # corresponding risk group factor
  
  # Call the plotting function
  plot_family_group(pc1, pc2, risk_vec, fam, order_name)
}

cat("Morphospace plots for all families have been generated!\n")