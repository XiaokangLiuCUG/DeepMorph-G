# Created by Xiaokang Liu on 31/10/2026 => xkliu@cug.edu.cn

library(data.table)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(parallel)
library(grid)

occ_path   <- "../Torres-2024-cleaneddata/6c086a1c-71b3-46a1-86e6-4d528f35a429/Dataset_Global-gastro_final_georeferenced_zscore3.csv"
morph_path <- "../ModernGastropods/PCA_data.rds"
tax_path   <- "../gastropods_taxonomy_combined.xlsx"
out_dir    <- "./ModernGastropods"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(occ_path))   stop("Occurrence file not found: ", occ_path)
if (!file.exists(morph_path)) stop("PCA_data.rds not found: ", morph_path)
if (!file.exists(tax_path))   stop("Taxonomy file not found: ", tax_path)

# ---------------------------- 1. Morphology & taxonomy ------------------------
pca_obj <- readRDS(morph_path)
if (inherits(pca_obj, "prcomp") || inherits(pca_obj, "gm.prcomp")) {
  morph <- pca_obj$x
} else {
  morph <- as.matrix(pca_obj)
}
if (is.null(rownames(morph))) stop("Morph matrix rownames are NULL.")
rownames(morph) <- trimws(tolower(rownames(morph)))
cat("morph species count:", nrow(morph), "\n")

taxonomy <- read_excel(tax_path)
taxonomy$species_clean <- trimws(tolower(taxonomy$species_clean))
taxonomy <- taxonomy[!duplicated(taxonomy$species_clean), ]
cat("taxonomy species count (deduplicated):", nrow(taxonomy), "\n")

calc_metrics <- function(sp_vec, mat, pc_use = 1:10) {
  idx <- which(rownames(mat) %in% sp_vec)
  n <- length(idx)
  if (n < 2) return(c(sum_var = NA_real_, sum_range = NA_real_))
  mat_use <- mat[idx, pc_use[pc_use <= ncol(mat)], drop = FALSE]
  if (ncol(mat_use) == 0) return(c(sum_var = NA_real_, sum_range = NA_real_))
  sum_var   <- sum(apply(mat_use, 2, var, na.rm = TRUE), na.rm = TRUE)
  sum_range <- sum(apply(mat_use, 2, function(x) diff(range(x, na.rm = TRUE))),
                   na.rm = TRUE)
  c(sum_var = sum_var, sum_range = sum_range)
}


# Grid level
occ <- fread(occ_path, select = c("Species", "Latitude", "Longitude"))
setnames(occ, c("species", "latitude", "longitude"))
occ <- occ[complete.cases(occ)]
occ <- occ[latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180]
occ[, species := trimws(tolower(species))]
cat("occ record count:", nrow(occ), "\n")

common_mo <- intersect(rownames(morph), occ$species)
cat("Grid-level species pool (morph ∩ occ):", length(common_mo), "\n")
if (length(common_mo) < 3) stop("Too few species in morph ∩ occ.")
occ_mo <- occ[species %in% common_mo]

occ_sf <- st_as_sf(occ_mo, coords = c("longitude", "latitude"), crs = 4326)
cell_size <- 5
grid <- st_make_grid(occ_sf, cellsize = cell_size, square = FALSE, flat_topped = TRUE)
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

sf_use_s2(FALSE)
valid <- st_is_valid(grid_sf)
if (!all(valid)) {
  grid_sf <- grid_sf[valid, ]
  grid_sf$grid_id <- 1:nrow(grid_sf)
}

intersections <- st_intersects(occ_sf, grid_sf)
occ_mo$cell_id <- sapply(intersections, function(x) if (length(x) > 0) x[1] else NA)
occ_mo <- occ_mo[!is.na(cell_id)]

grid_summary <- occ_mo[, .(
  occurrence_count = .N,
  species_list     = list(unique(species))
), by = cell_id]
setnames(grid_summary, "cell_id", "grid_id")
grid_summary[, richness := lengths(species_list)]

grid_sf <- merge(grid_sf, grid_summary, by = "grid_id", all.x = TRUE)
grid_sf <- grid_sf[!is.na(grid_sf$richness) & grid_sf$richness > 10, ]
cat("Grid-level retained grids (richness > 10):", nrow(grid_sf), "\n")

sp_lists <- grid_sf$species_list
ncores <- max(1, detectCores() - 1)

if (Sys.info()["sysname"] == "Windows") {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("calc_metrics", "morph"), envir = environment())
  res <- parLapply(cl, sp_lists, function(sp) calc_metrics(sp, morph))
  stopCluster(cl)
} else {
  res <- mclapply(sp_lists, function(sp) calc_metrics(sp, morph), mc.cores = ncores)
}
grid_sf$sum_var   <- sapply(res, function(x) x["sum_var"])
grid_sf$sum_range <- sapply(res, function(x) x["sum_range"])

grid_df <- st_set_geometry(grid_sf, NULL) %>%
  filter(!is.na(sum_var) & !is.na(sum_range)) %>%
  mutate(log_richness = log10(richness))


# Family level (morph ∩ taxonomy only, no occ needed)

common_tax <- intersect(rownames(morph), taxonomy$species_clean)
cat("Family-level species pool (morph ∩ taxonomy):", length(common_tax), "\n")
if (length(common_tax) < 3) stop("Too few species in morph ∩ taxonomy.")

morph_tax <- morph[common_tax, , drop = FALSE]
tax_sub   <- taxonomy[taxonomy$species_clean %in% common_tax, ]
tax_sub   <- tax_sub[!is.na(tax_sub$Family) & tax_sub$Family != "", ]
cat("Family-level usable species (non-empty Family):", nrow(tax_sub), "\n")

family_list <- split(tax_sub$species_clean, tax_sub$Family)
family_list <- lapply(family_list, function(sp) intersect(sp, common_tax))
family_list <- family_list[lengths(family_list) > 0]

family_df_all <- do.call(rbind, lapply(names(family_list), function(fam) {
  sp <- family_list[[fam]]
  m  <- calc_metrics(sp, morph_tax)
  data.frame(family    = fam,
             richness  = length(sp),
             sum_var   = as.numeric(m["sum_var"]),
             sum_range = as.numeric(m["sum_range"]),
             stringsAsFactors = FALSE)
}))
cat("Number of families with richness >= 2:", nrow(family_df_all), "\n")

family_df <- family_df_all %>%
  filter(!is.na(sum_var) & !is.na(sum_range) & richness > 10) %>%
  mutate(log_richness = log10(richness))
cat("Family-level retained families (richness > 10):", nrow(family_df), "\n")
cat("Family-level richness distribution:\n"); print(summary(family_df$richness))
cat("Top 10 families by richness:\n")
print(head(family_df[order(-family_df$richness),
                     c("family", "richness", "sum_var", "sum_range")], 10))

x_breaks <- log10(c(10, 100, 1000, 10000))
x_label_fun <- function(x) {
  format(round(10^x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fit_linear <- function(df, x_col, y_col, x_label, y_label, title) {
  df_fit <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), ]
  fit <- lm(df_fit[[y_col]] ~ df_fit[[x_col]])
  s   <- summary(fit)
  r2  <- s$r.squared
  r   <- cor(df_fit[[x_col]], df_fit[[y_col]])
  p   <- s$coefficients[2, 4]
  label <- sprintf("R\u00b2 = %.3f\nr = %.3f\np = %s\nn = %d",
                   r2, r, format.pval(p, digits = 3), nrow(df_fit))
  
  p_plot <- ggplot(df_fit, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red",
                fill = "pink", alpha = 0.3) +
    annotate("text", x = Inf, y = Inf, label = label,
             hjust = 1.1, vjust = 1.1, size = 4.2) +
    scale_x_continuous(breaks = x_breaks, labels = x_label_fun) +
    labs(x = x_label, y = y_label, title = title) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          panel.grid.minor = element_blank())
  
  list(plot = p_plot, r2 = r2, r = r, p = p, n = nrow(df_fit), fit = fit)
}

fg_range <- fit_linear(grid_df, "log_richness", "sum_range",
                       "Richness", "Sum of ranges (PC1-10)",
                       "Grid level: Richness vs Sum_range")
fg_var   <- fit_linear(grid_df, "log_richness", "sum_var",
                       "Richness", "Sum of variances (PC1-10)",
                       "Grid level: Richness vs Sum_var")
ff_range <- fit_linear(family_df, "log_richness", "sum_range",
                       "Richness (per family)", "Sum of ranges (PC1-10)",
                       "Family level: Richness vs Sum_range")
ff_var   <- fit_linear(family_df, "log_richness", "sum_var",
                       "Richness (per family)", "Sum of variances (PC1-10)",
                       "Family level: Richness vs Sum_var")

results_tab <- data.frame(
  level    = c("grid",   "grid",   "family", "family"),
  relation = c("log10(Richness)~Sum_range", "log10(Richness)~Sum_var",
               "log10(Richness)~Sum_range", "log10(Richness)~Sum_var"),
  n        = c(fg_range$n, fg_var$n, ff_range$n, ff_var$n),
  R2       = c(fg_range$r2, fg_var$r2, ff_range$r2, ff_var$r2),
  r        = c(fg_range$r,  fg_var$r,  ff_range$r,  ff_var$r),
  p        = c(fg_range$p,  fg_var$p,  ff_range$p,  ff_var$p)
)
print(results_tab, row.names = FALSE)
write.csv(results_tab,
          file.path(out_dir, "linear_fit_logRichness_range_var.csv"),
          row.names = FALSE)

fwrite(as.data.table(grid_df),
       file.path(out_dir, "grid_level_data.csv"))
fwrite(as.data.table(family_df),
       file.path(out_dir, "family_level_data.csv"))
fwrite(as.data.table(family_df_all),
       file.path(out_dir, "family_level_all_families.csv"))

pdf(file.path(out_dir, "LinearFit_logRichness_range_var.pdf"),
    width = 12, height = 10)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 2)))
print(fg_range$plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(fg_var$plot,   vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ff_range$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(ff_var$plot,   vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 2)))
print(fg_range$plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(fg_var$plot,   vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ff_range$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(ff_var$plot,   vp = viewport(layout.pos.row = 2, layout.pos.col = 2))