# ==============================================================================
# Biocontrol Suitability + Management Gap Mapping for India
# Low-threshold early-warning version preserving 4-class overlap structure
# Updated: 2026-08-05
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(maxnet)
  library(tidyterra)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
  library(grid)
})

# -----------------------------------------------------------------------------
# 0. Paths and directories
# -----------------------------------------------------------------------------
base_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"
bio_out_dir <- file.path(base_dir, "biocontrol")
fig_out <- file.path(base_dir, "biocontrol_4class_forecast_lowthreshold.png")
summary_out <- file.path(base_dir, "biocontrol_4class_area_summary_lowthreshold.csv")
curve_out <- file.path(base_dir, "zygogramma_thermal_curve.png")
parth_only_curr_out <- file.path(base_dir, "parthenium_only_current_lowthreshold.tif")
parth_only_2030_out <- file.path(base_dir, "parthenium_only_2030_lowthreshold.tif")
parth_only_2050_out <- file.path(base_dir, "parthenium_only_2050_lowthreshold.tif")
parth_only_points_2030 <- file.path(base_dir, "parthenium_only_2030_points_lowthreshold.gpkg")
parth_only_points_2050 <- file.path(base_dir, "parthenium_only_2050_points_lowthreshold.gpkg")

if (!dir.exists(bio_out_dir)) dir.create(bio_out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Zygogramma thermal suitability function
# -----------------------------------------------------------------------------
# Broad climatic permissiveness proxy based on annual mean temperature.
# This is not a mechanistic establishment model.
zygogramma_suitability <- function(temp, Tmin = 10, Topt = 27, Tmax = 38) {
  if (!(Tmin < Topt && Topt < Tmax)) {
    stop("Thermal parameters must satisfy Tmin < Topt < Tmax")
  }
  
  exponent <- (Tmax - Topt) / (Topt - Tmin)
  
  suit <- ((temp - Tmin) * (Tmax - temp)^exponent) /
    ((Topt - Tmin) * (Tmax - Topt)^exponent)
  
  suit[temp <= Tmin | temp >= Tmax] <- 0
  suit[!is.finite(suit)] <- NA
  suit[suit < 0] <- 0
  suit[suit > 1] <- 1
  return(suit)
}

# -----------------------------------------------------------------------------
# 2. Diagnostic thermal curve
# -----------------------------------------------------------------------------
temp_check <- seq(0, 45, by = 0.1)
suit_check <- zygogramma_suitability(temp_check)

png(curve_out, width = 800, height = 600, res = 120)
plot(
  temp_check, suit_check,
  type = "l", lwd = 2.5, col = "forestgreen",
  xlab = "Mean Temperature (°C)",
  ylab = "Thermal Suitability Index (0-1)",
  main = expression(paste("Thermal Suitability Curve for ", italic("Zygogramma bicolorata")))
)
abline(v = c(10, 38), lty = 3, col = "firebrick", lwd = 1.5)
abline(v = 27, lty = 2, col = "blue", lwd = 1.5)
legend(
  "topright",
  legend = c("Suitability Index", "Topt (27°C)", "Tmin/Tmax (10°C / 38°C)"),
  col = c("forestgreen", "blue", "firebrick"),
  lty = c(1, 2, 3), lwd = c(2.5, 1.5, 1.5), bty = "n"
)
dev.off()

# -----------------------------------------------------------------------------
# 3. Load climate rasters
# -----------------------------------------------------------------------------
mat_curr <- rast(file.path(base_dir, "climate/processed/mean_ann_temp.tif"))
mat_2030_all <- rast(file.path(base_dir, "climate/future/env_2030.tif"))
mat_2050_all <- rast(file.path(base_dir, "climate/future/env_2050.tif"))

if (!("mean_ann_temp" %in% names(mat_2030_all))) {
  stop("Layer 'mean_ann_temp' not found in env_2030.tif")
}
if (!("mean_ann_temp" %in% names(mat_2050_all))) {
  stop("Layer 'mean_ann_temp' not found in env_2050.tif")
}

mat_2030 <- mat_2030_all[["mean_ann_temp"]]
mat_2050 <- mat_2050_all[["mean_ann_temp"]]

mat_2030 <- project(mat_2030, mat_curr, method = "bilinear")
mat_2050 <- project(mat_2050, mat_curr, method = "bilinear")

# -----------------------------------------------------------------------------
# 4. Generate Zygogramma climatic suitability rasters
# -----------------------------------------------------------------------------
bio_suit_curr <- zygogramma_suitability(mat_curr)
bio_suit_2030 <- zygogramma_suitability(mat_2030)
bio_suit_2050 <- zygogramma_suitability(mat_2050)

writeRaster(bio_suit_curr, file.path(bio_out_dir, "bio_suit_curr.tif"), overwrite = TRUE)
writeRaster(bio_suit_2030, file.path(bio_out_dir, "bio_suit_2030.tif"), overwrite = TRUE)
writeRaster(bio_suit_2050, file.path(bio_out_dir, "bio_suit_2050.tif"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 5. Load Maxnet model and predictor stacks
# -----------------------------------------------------------------------------
maxnet_predict <- function(model, data, ...) {
  predict(model, newdata = data, type = "cloglog", ...)
}

best_mod <- readRDS(file.path(base_dir, "best_maxnet_model.rds"))

envs_curr <- rast(c(
  file.path(base_dir, "climate/processed/mean_ann_temp.tif"),
  file.path(base_dir, "climate/processed/total_ann_precip.tif"),
  file.path(base_dir, "climate/processed/temp_seasonality.tif")
))

env_2030 <- rast(file.path(base_dir, "climate/future/env_2030.tif"))
env_2050 <- rast(file.path(base_dir, "climate/future/env_2050.tif"))

required_vars <- names(envs_curr)
if (!all(required_vars %in% names(env_2030))) {
  stop("env_2030 is missing one or more required predictors: ", paste(required_vars, collapse = ", "))
}
if (!all(required_vars %in% names(env_2050))) {
  stop("env_2050 is missing one or more required predictors: ", paste(required_vars, collapse = ", "))
}

env_2030 <- env_2030[[required_vars]]
env_2050 <- env_2050[[required_vars]]

env_2030 <- project(env_2030, envs_curr, method = "bilinear")
env_2050 <- project(env_2050, envs_curr, method = "bilinear")

# -----------------------------------------------------------------------------
# 6. Predict Parthenium suitability
# -----------------------------------------------------------------------------
cat("Predicting Parthenium suitability across climate scenarios...\n")
part_curr <- predict(envs_curr, best_mod, fun = maxnet_predict, na.rm = TRUE)
part_2030 <- predict(env_2030, best_mod, fun = maxnet_predict, na.rm = TRUE)
part_2050 <- predict(env_2050, best_mod, fun = maxnet_predict, na.rm = TRUE)

writeRaster(part_curr, file.path(bio_out_dir, "parthenium_suit_curr.tif"), overwrite = TRUE)
writeRaster(part_2030, file.path(bio_out_dir, "parthenium_suit_2030.tif"), overwrite = TRUE)
writeRaster(part_2050, file.path(bio_out_dir, "parthenium_suit_2050.tif"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 7. Thresholds
# -----------------------------------------------------------------------------
# Early-warning / invasion-pathway framing:
# - Use a lower Parthenium threshold derived from a permissive occurrence-based
#   cutoff to reveal marginally suitable expansion space.
# - Keep the insect threshold stricter because it is only a thermal proxy.
p_thresh <- 0.4952
bio_thresh <- 0.75

part_bin_curr <- ifel(part_curr >= p_thresh, 1, 0)
part_bin_2030 <- ifel(part_2030 >= p_thresh, 1, 0)
part_bin_2050 <- ifel(part_2050 >= p_thresh, 1, 0)

bio_bin_curr <- ifel(bio_suit_curr >= bio_thresh, 1, 0)
bio_bin_2030 <- ifel(bio_suit_2030 >= bio_thresh, 1, 0)
bio_bin_2050 <- ifel(bio_suit_2050 >= bio_thresh, 1, 0)

writeRaster(part_bin_curr, file.path(bio_out_dir, "part_bin_curr_lowthreshold.tif"), overwrite = TRUE)
writeRaster(part_bin_2030, file.path(bio_out_dir, "part_bin_2030_lowthreshold.tif"), overwrite = TRUE)
writeRaster(part_bin_2050, file.path(bio_out_dir, "part_bin_2050_lowthreshold.tif"), overwrite = TRUE)

writeRaster(bio_bin_curr, file.path(bio_out_dir, "bio_bin_curr_strict.tif"), overwrite = TRUE)
writeRaster(bio_bin_2030, file.path(bio_out_dir, "bio_bin_2030_strict.tif"), overwrite = TRUE)
writeRaster(bio_bin_2050, file.path(bio_out_dir, "bio_bin_2050_strict.tif"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 8. Common masks
# -----------------------------------------------------------------------------
common_mask_curr <- !is.na(part_curr) & !is.na(bio_suit_curr)
common_mask_2030 <- !is.na(part_2030) & !is.na(bio_suit_2030)
common_mask_2050 <- !is.na(part_2050) & !is.na(bio_suit_2050)

part_bin_curr <- mask(part_bin_curr, common_mask_curr, maskvalues = 0)
part_bin_2030 <- mask(part_bin_2030, common_mask_2030, maskvalues = 0)
part_bin_2050 <- mask(part_bin_2050, common_mask_2050, maskvalues = 0)

bio_bin_curr <- mask(bio_bin_curr, common_mask_curr, maskvalues = 0)
bio_bin_2030 <- mask(bio_bin_2030, common_mask_2030, maskvalues = 0)
bio_bin_2050 <- mask(bio_bin_2050, common_mask_2050, maskvalues = 0)

# -----------------------------------------------------------------------------
# 9. Four-class map creation
# -----------------------------------------------------------------------------
create_4class_code_map <- function(part_bin, bio_bin, mask_template) {
  code_rast <- classify(
    part_bin * 2 + bio_bin,
    rcl = matrix(c(
      0, 1,
      1, 4,
      2, 3,
      3, 2
    ), ncol = 2, byrow = TRUE)
  )
  code_rast <- mask(code_rast, mask_template, maskvalues = 0)
  return(code_rast)
}

make_factor_map <- function(code_rast) {
  risk_rast <- as.factor(code_rast)
  levels(risk_rast) <- data.frame(
    value = 1:4,
    category = c(
      "Unsuitable for both",
      "Overlap: both suitable",
      "Parthenium only (management gap)",
      "Biocontrol proxy only"
    )
  )
  return(risk_rast)
}

risk_code_curr <- create_4class_code_map(part_bin_curr, bio_bin_curr, common_mask_curr)
risk_code_2030 <- create_4class_code_map(part_bin_2030, bio_bin_2030, common_mask_2030)
risk_code_2050 <- create_4class_code_map(part_bin_2050, bio_bin_2050, common_mask_2050)

risk_curr <- make_factor_map(risk_code_curr)
risk_2030 <- make_factor_map(risk_code_2030)
risk_2050 <- make_factor_map(risk_code_2050)

writeRaster(risk_code_curr, file.path(bio_out_dir, "risk_4class_curr_lowthreshold.tif"), overwrite = TRUE)
writeRaster(risk_code_2030, file.path(bio_out_dir, "risk_4class_2030_lowthreshold.tif"), overwrite = TRUE)
writeRaster(risk_code_2050, file.path(bio_out_dir, "risk_4class_2050_lowthreshold.tif"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 10. Explicit Parthenium-only rasters and point exports
# -----------------------------------------------------------------------------
parth_only_curr <- ifel(part_bin_curr == 1 & bio_bin_curr == 0, 1, NA)
parth_only_2030 <- ifel(part_bin_2030 == 1 & bio_bin_2030 == 0, 1, NA)
parth_only_2050 <- ifel(part_bin_2050 == 1 & bio_bin_2050 == 0, 1, NA)

writeRaster(parth_only_curr, parth_only_curr_out, overwrite = TRUE)
writeRaster(parth_only_2030, parth_only_2030_out, overwrite = TRUE)
writeRaster(parth_only_2050, parth_only_2050_out, overwrite = TRUE)

parth_only_2030_pts <- as.points(parth_only_2030, values = TRUE, na.rm = TRUE)
if (!is.null(parth_only_2030_pts) && nrow(parth_only_2030_pts) > 0) {
  writeVector(parth_only_2030_pts, parth_only_points_2030, overwrite = TRUE)
}

parth_only_2050_pts <- as.points(parth_only_2050, values = TRUE, na.rm = TRUE)
if (!is.null(parth_only_2050_pts) && nrow(parth_only_2050_pts) > 0) {
  writeVector(parth_only_2050_pts, parth_only_points_2050, overwrite = TRUE)
}

# -----------------------------------------------------------------------------
# 11. Area summary using integer-coded raster
# -----------------------------------------------------------------------------
calc_area_stats <- function(risk_code_rast, time_label) {
  area_r <- cellSize(risk_code_rast, unit = "km")
  z <- zonal(area_r, risk_code_rast, fun = "sum", na.rm = TRUE)
  colnames(z) <- c("value", "area_km2")
  
  full_df <- data.frame(
    value = 1:4,
    category = c(
      "Unsuitable for both",
      "Overlap: both suitable",
      "Parthenium only (management gap)",
      "Biocontrol proxy only"
    )
  )
  
  stat_df <- merge(full_df, z, by = "value", all.x = TRUE)
  stat_df$area_km2[is.na(stat_df$area_km2)] <- 0
  
  total_area <- sum(stat_df$area_km2, na.rm = TRUE)
  stat_df$pct <- if (total_area > 0) (stat_df$area_km2 / total_area) * 100 else 0
  stat_df$Period <- time_label
  stat_df <- stat_df[, c("Period", "value", "category", "area_km2", "pct")]
  return(stat_df)
}

stats_curr <- calc_area_stats(risk_code_curr, "Current baseline")
stats_2030 <- calc_area_stats(risk_code_2030, "2030 (SSP2-4.5)")
stats_2050 <- calc_area_stats(risk_code_2050, "2050 (SSP2-4.5)")

summary_all <- rbind(stats_curr, stats_2030, stats_2050)
write.csv(summary_all, summary_out, row.names = FALSE)
print(summary_all)

# -----------------------------------------------------------------------------
# 12. Colorblind-safe palette and plotting
# -----------------------------------------------------------------------------
risk_cols <- c(
  "Unsuitable for both" = "#BDBDBD",
  "Overlap: both suitable" = "#0072B2",
  "Parthenium only (management gap)" = "#D55E00",
  "Biocontrol proxy only" = "#009E73"
)

plot_risk_panel <- function(risk_rast, title_str, add_annotations = FALSE) {
  p <- ggplot() +
    geom_spatraster(data = risk_rast) +
    scale_fill_manual(
      values = risk_cols,
      na.translate = FALSE,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 5))
    ) +
    theme_bw() +
    labs(title = title_str) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA),
      legend.text = element_text(size = 10),
      legend.spacing.x = unit(0.15, "cm")
    )
  
  if (add_annotations) {
    p <- p +
      annotation_scale(location = "bl", width_hint = 0.25, style = "ticks", text_cex = 0.7) +
      annotation_north_arrow(
        location = "tl",
        which_north = "true",
        pad_x = unit(0.2, "in"),
        pad_y = unit(0.2, "in"),
        style = north_arrow_minimal(text_size = 8)
      )
  }
  
  return(p)
}

p1 <- plot_risk_panel(risk_curr, "Current baseline", add_annotations = TRUE)
p2 <- plot_risk_panel(risk_2030, "2030 projection (SSP2-4.5)")
p3 <- plot_risk_panel(risk_2050, "2050 projection (SSP2-4.5)")

final_figure <- (p1 + p2 + p3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(2, 2, 2, 2),
    legend.box.margin = margin(0, 0, 0, 0)
  )

final_figure <- final_figure +
  plot_annotation(
    title = "Biocontrol climatic suitability vs. Parthenium suitability forecasts — India",
    subtitle = sprintf(
      "Parthenium threshold = %.4f | Zygogramma climatic suitability threshold = %.2f",
      p_thresh, bio_thresh
    ),
    caption = paste(
      "Blue = overlap; Red-orange = Parthenium-only management gap;",
      "Green = biocontrol proxy only; Grey = unsuitable for both.",
      "Interpret insect suitability as a broad thermal proxy, not proof of establishment."
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0),
      plot.caption = element_text(size = 9.5, color = "grey35", hjust = 0)
    )
  )

ggsave(
  filename = fig_out,
  plot = final_figure,
  width = 16,
  height = 7.2,
  dpi = 300
)

cat("Biocontrol 4-class low-threshold map and summary table saved successfully.\n")
