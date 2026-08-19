# ==============================================================================
# SCRIPT 02: CANDIDATE GAM MODEL FITTING & RDS EXPORT
# ==============================================================================

library(tidyverse)
library(mgcv)
library(sf)

# 1. Load survey data and factorize variables
survey_raw <- read.csv("/home/anaga-ambady/Documents/DISSERT/DATA/MODEL/results/survey_data.csv") %>%
  mutate(
    hab_group = as.factor(hab_group),
    road_type = as.factor(road_type),
    plotID    = as.factor(plotID)
  )

# 2. Reproject coordinates to UTM Zone 42N (EPSG:32642) and scale to kilometers
survey_sf <- survey_raw %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 32642)

coords_m <- st_coordinates(survey_sf)

survey_data <- survey_raw %>%
  mutate(
    x_km = coords_m[, 1] / 1000,
    y_km = coords_m[, 2] / 1000
  )

rds_dir <- "/home/anaga-ambady/Documents/DISSERT/EDITS/section 1/gam_rds"
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# --- M5: Full Primary Model (REML, k=200 GP spatial smooth) ---
gam_gdd_k_200 <- gam(
  presence ~ hab_group + road_type +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, by = hab_group, k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 200),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M4: Spatial Reduced Model (REML, k=100 GP spatial smooth) ---
gam_gdd_k_100 <- gam(
  presence ~ hab_group + road_type +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, by = hab_group, k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 100),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M3: Altitude-Removed Model (REML, k=100) ---
gam_gdd_xalt <- gam(
  presence ~ hab_group + road_type +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, by = hab_group, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 100),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M2: Factor-Smooth Model (REML, k=100) ---
gam_gdd_fs <- gam(
  presence ~ hab_group + road_type +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, hab_group, bs = "fs", k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 100),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M1: Non-Spatial Model (REML) ---
gam_gdd_no_loc <- gam(
  presence ~ hab_group + road_type +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, by = hab_group, k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M7: No Road-Type Model (REML) ---
gam_gdd_noroad <- gam(
  presence ~ hab_group +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, by = hab_group, k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 100),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# --- M6: No Road-Type Factor-Smooth Model (REML) ---
gam_gdd_noroad_sm <- gam(
  presence ~ hab_group +
    s(precip_sum_30d, k = 25) +
    s(gdd_30d, k = 25) +
    s(gdd_30d, hab_group, bs = "fs", k = 15) +
    s(altitude, k = 15) +
    s(hab_perc, k = 10) +
    s(x_km, y_km, bs = "gp", k = 100),
  family = binomial(link = "logit"), method = "REML", data = survey_data
)

# Save RDS Objects
saveRDS(gam_gdd_k_200,      file.path(rds_dir, "gam_gdd_k_200.rds"))
saveRDS(gam_gdd_k_100,      file.path(rds_dir, "gam_gdd_k_100.rds"))
saveRDS(gam_gdd_xalt,       file.path(rds_dir, "gam_gdd_xalt.rds"))
saveRDS(gam_gdd_fs,         file.path(rds_dir, "gam_gdd_fs.rds"))
saveRDS(gam_gdd_no_loc,     file.path(rds_dir, "gam_gdd_no_loc.rds"))
saveRDS(gam_gdd_noroad,     file.path(rds_dir, "gam_gdd_noroad.rds"))
saveRDS(gam_gdd_noroad_sm,  file.path(rds_dir, "gam_gdd_noroad_sm.rds"))

graphics.off()           # Close open graphic devices
par(mar = c(3, 3, 2, 1)) # Set smaller inner margins (bottom, left, top, right)
gam.check(gam_gdd_k_200)