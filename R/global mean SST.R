# Created by Xiaokang Liu on 15/10/2026 => xkliu@cug.edu.cn
library(ncdf4)
library(tidyverse)
library(lubridate)

# ===== 3. Read data =====
nc_file <- "sst.mnmean.nc"
nc <- nc_open(nc_file)

sst_var  <- ncvar_get(nc, "sst")
lat      <- ncvar_get(nc, "lat")
lon      <- ncvar_get(nc, "lon")
time_raw <- ncvar_get(nc, "time")

time_units  <- ncatt_get(nc, "time", "units")$value
origin_str  <- sub(".*since ", "", time_units)
time_origin <- as.Date(origin_str, format = "%Y-%m-%d")
time        <- as.Date(time_raw, origin = time_origin)

nc_close(nc)

# ===== 4. Bootstrap zonal mean =====
# Resampling unit: years (with replacement)
# Returns: mean, lower (2.5%), upper (97.5%)
boot_zonal_mean <- function(start_year, end_year,
                            n_boot = 1000,
                            probs  = c(0.025, 0.975),
                            seed   = 123) {
  
  set.seed(seed)  # reproducible
  
  year_series <- year(time)
  idx <- which(year_series >= start_year & year_series <= end_year)
  if (length(idx) == 0) stop("No data found for this period")
  
  sst_sub   <- sst_var[, , idx, drop = FALSE]
  time_sub  <- time[idx]
  years_sub <- year(time_sub)
  unique_years <- sort(unique(years_sub))
  n_year <- length(unique_years)
  
  # --- First compute the "zonal mean field" for each year ---
  # Result dimensions: [n_lat, n_year]
  zonal_yearly <- matrix(NA, nrow = length(lat), ncol = n_year)
  for (i in seq_len(n_year)) {
    y <- unique_years[i]
    idx_y <- which(years_sub == y)
    
    # Annual mean field for that year [lon, lat]
    sst_y <- apply(sst_sub[, , idx_y, drop = FALSE],
                   MARGIN = c(1, 2), mean, na.rm = TRUE)
    
    # Average along longitude -> zonal mean [lat]
    zonal_yearly[, i] <- apply(sst_y, 2, mean, na.rm = TRUE)
  }
  
  # --- Bootstrap: resample years with replacement n_year times, repeat n_boot times ---
  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(lat))
  for (b in seq_len(n_boot)) {
    samp <- sample(seq_len(n_year), size = n_year, replace = TRUE)
    boot_mat[b, ] <- apply(zonal_yearly[, samp, drop = FALSE],
                           1, mean, na.rm = TRUE)
  }
  
  # --- Mean & quantile interval ---
  zonal_mean <- apply(zonal_yearly, 1, mean, na.rm = TRUE)
  ci <- apply(boot_mat, 2, quantile, probs = probs, na.rm = TRUE)
  
  data.frame(
    latitude = lat,
    sst      = zonal_mean,
    lower    = ci[1, ],
    upper    = ci[2, ]
  )
}

# ===== 5. Two periods =====
df1 <- boot_zonal_mean(1854, 1860, n_boot = 1000)
df2 <- boot_zonal_mean(1971, 2000, n_boot = 1000)

df1$period <- "1854-1860"
df2$period <- "1971-2000"
df <- rbind(df1, df2)
df$period <- factor(df$period, levels = c("1854-1860", "1971-2000"))

# ===== 6. Plot =====
p <- ggplot(df, aes(x = latitude, y = sst, color = period, fill = period)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("1854-1860" = "darkred",
                                "1971-2000" = "steelblue"), name = "Period") +
  scale_fill_manual(values  = c("1854-1860" = "darkred",
                                "1971-2000" = "steelblue"), name = "Period") +
  labs(
    title = "Zonal Mean SST with Bootstrap 95% CI (1000 resamples)",
    x = "Latitude (°N)",
    y = "SST (°C)"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, 30)) +
  theme(legend.position = "bottom")

print(p)

ggsave("zonal_sst_bootstrap_CI.pdf", plot = p,
       width = 10, height = 6, units = "in")