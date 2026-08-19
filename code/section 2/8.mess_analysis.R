# ==============================================================================
# Multivariate Environmental Similarity Surface (MESS) Analysis
# Script 8 in the India Modeling Pipeline
#
# Flags extrapolation risk in future climate projections (2030, 2050)
# relative to the calibration training space (presence + background points).
# Updated to align with current pipeline paths and plotting standards.
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(predicts)
  library(tidyterra)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
})

# -----------------------------------------------------------------------------
# 1. Paths and Directories
# -----------------------------------------------------------------------------
base_dir      <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"
processed_dir <- file.path(base_dir, "climate/processed")
future_dir    <- file.path(base_dir, "climate/future")

fig_out       <- file.path(future_dir, "mess_analysis_future.png")
mess_2030_out <- file.path(future_dir, "mess_2030.tif")
mess_2050_out <- file.path(future_dir, "mess_2050.tif")

# -----------------------------------------------------------------------------
# 2. Load Calibration Rasters and Point Data
# -----------------------------------------------------------------------------
envs_curr <- c(
  rast(file.path(processed_dir, "mean_ann_temp.tif")),
  rast(file.path(processed_dir, "total_ann_precip.tif")),
  rast(file.path(processed_dir, "temp_seasonality.tif"))
)
names(envs_curr) <- c("mean_ann_temp", "total_ann_precip", "temp_seasonality")

occs_clean <- read.csv(file.path(base_dir, "parthenium_thinned_clean.csv"))
bg_points  <- read.csv(file.path(base_dir, "background_points.csv"))

# Determine coordinate column names dynamically
get_coords <- function(df) {
  if (all(c("lon", "lat") %in% names(df))) return(df[, c("lon", "lat")])
  if (all(c("longitude", "latitude") %in% names(df))) return(df[, c("longitude", "latitude")])
  if (all(c("x", "y") %in% names(df))) return(df[, c("x", "y")])
  stop("Coordinate columns not recognized in input data frame.")
}

occs_coords <- get_coords(occs_clean)
bg_coords   <- get_coords(bg_points)

# Extract environmental values across calibration domain
occs_vals <- terra::extract(envs_curr, occs_coords)[, -1]
bg_vals   <- terra::extract(envs_curr, bg_coords)[, -1]

training_envelope <- na.omit(rbind(occs_vals, bg_vals))

# -----------------------------------------------------------------------------
# 3. Load and Align Future Climate Rasters
# -----------------------------------------------------------------------------
env_2030_raw <- rast(file.path(future_dir, "env_2030.tif"))
env_2050_raw <- rast(file.path(future_dir, "env_2050.tif"))

required_vars <- names(envs_curr)
if (!all(required_vars %in% names(env_2030_raw))) {
  stop("env_2030 is missing required predictors: ", paste(setdiff(required_vars, names(env_2030_raw)), collapse = ", "))
}
if (!all(required_vars %in% names(env_2050_raw))) {
  stop("env_2050 is missing required predictors: ", paste(setdiff(required_vars, names(env_2050_raw)), collapse = ", "))
}

env_2030 <- env_2030_raw[[required_vars]]
env_2050 <- env_2050_raw[[required_vars]]

env_2030 <- project(env_2030, envs_curr, method = "bilinear")
env_2050 <- project(env_2050, envs_curr, method = "bilinear")

# -----------------------------------------------------------------------------
# 4. Compute MESS Rasters & Spatial Masking
# -----------------------------------------------------------------------------
cat("Calculating MESS rasters for 2030 and 2050 scenarios...\n")
mess_2030 <- predicts::mess(env_2030, training_envelope)
mess_2050 <- predicts::mess(env_2050, training_envelope)

# Mask using the baseline raster template to eliminate floating eastern artifacts
mess_2030 <- mask(mess_2030, envs_curr[[1]])
mess_2050 <- mask(mess_2050, envs_curr[[1]])

writeRaster(mess_2030, mess_2030_out, overwrite = TRUE)
writeRaster(mess_2050, mess_2050_out, overwrite = TRUE)

# Calculate global minimum/maximum values across both projections for a unified scale
global_min <- min(c(minmax(mess_2030)[1], minmax(mess_2050)[1]), na.rm = TRUE)
global_max <- max(c(minmax(mess_2030)[2], minmax(mess_2050)[2]), na.rm = TRUE)

# Dynamic caption generation based on extrapolation detection
caption_text <- if (global_min < 0) {
  "Negative values (orange) indicate novel climatic conditions outside the training envelope."
} else {
  "Values >= 0 indicate future climate conditions remain entirely within the calibration training space."
}

# -----------------------------------------------------------------------------
# 5. Publication-Quality Visualization
# -----------------------------------------------------------------------------
plot_mess_panel <- function(mess_rast, title_str, limits_val, add_annotations = FALSE) {
  p <- ggplot() +
    geom_spatraster(data = mess_rast) +
    scale_fill_gradient2(
      low = "#D55E00",      # Negative values: Novel climate space (extrapolation risk)
      mid = "#F0E442",      # Near zero: Marginal conditions
      high = "#0072B2",     # Positive values: Well within training envelope
      midpoint = 0,
      limits = limits_val,  # Forces identical color limits for legend unification
      na.value = "transparent",
      name = "MESS Index"
    ) +
    theme_bw() +
    labs(title = title_str) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold")
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

limits_range <- c(global_min, global_max)

p1 <- plot_mess_panel(mess_2030, "2030 Projection (SSP2-4.5)", limits_range, add_annotations = TRUE)
p2 <- plot_mess_panel(mess_2050, "2050 Projection (SSP2-4.5)", limits_range)

final_figure <- (p1 + p2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

final_figure <- final_figure +
  plot_annotation(
    title = "Multivariate Environmental Similarity Surface (MESS) Analysis",
    subtitle = "Extrapolation risk assessment relative to calibration training space",
    caption = caption_text,
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0),
      plot.caption = element_text(size = 9.5, color = "grey35", hjust = 0)
    )
  )

ggsave(
  filename = fig_out,
  plot = final_figure,
  width = 13,
  height = 6.5,
  dpi = 300
)

cat("MESS analysis complete. Plots and rasters successfully saved to climate/future/\n")