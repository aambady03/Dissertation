# =========================================================
# Step 3: Background (pseudo-absence) Point Sampling for MaxEnt
# Bridges bias surface output (calibration_area + bias_file_norm)
# to downstream model training steps.
# =========================================================

library(terra)

# -------------------------------------------------------------------------
# 1. Load inputs produced by calibration_zone.R and bias_file.R
# -------------------------------------------------------------------------
base_dir   <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"
calib_path <- file.path(base_dir, "calibration_area.shp")
bias_path  <- file.path(base_dir, "bias/results/bias_file.tif")

if (!file.exists(calib_path) || !file.exists(bias_path)) {
  stop("Missing calibration area shapefile or normalized bias raster. Run Step 1 and Step 2 first!")
}

calibration_area <- vect(calib_path)
bias_file_norm   <- rast(bias_path)

# Ensure output directory exists
out_dir <- file.path(base_dir, "bias/results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------------------------------------------------
# 2. Re-mask the bias surface to the calibration area
# -------------------------------------------------------------------------
# Check CRS alignment before masking
if (crs(bias_file_norm) != crs(calibration_area)) {
  calibration_area <- project(calibration_area, crs(bias_file_norm))
}

bias_masked <- mask(bias_file_norm, calibration_area)

# -------------------------------------------------------------------------
# 3. Convert bias raster to a data frame for sampling
# -------------------------------------------------------------------------
bias_df <- as.data.frame(bias_masked, xy = TRUE, na.rm = TRUE)
names(bias_df) <- c("lon", "lat", "weight")

# Filter out zero or negative weights to prevent sampling issues
bias_df <- bias_df[bias_df$weight > 0, ]

cat("Candidate background cells available:", nrow(bias_df), "\n")

# -------------------------------------------------------------------------
# 4. Weighted random sample of background points
# -------------------------------------------------------------------------
set.seed(42)  # Set seed for exact reproducibility
n_bg <- 10000

if (nrow(bias_df) < n_bg) {
  stop("Fewer candidate cells (", nrow(bias_df), ") than requested background points (",
       n_bg, "). Lower n_bg or check calibration area/resolution.")
}

bg_idx <- sample(
  x       = seq_len(nrow(bias_df)),
  size    = n_bg,
  replace = FALSE,
  prob    = bias_df$weight
)

bg_points <- bias_df[bg_idx, c("lon", "lat")]

cat("Background points drawn:", nrow(bg_points), "\n")

# -------------------------------------------------------------------------
# 5. Diagnostic plot — confirm background points track bias surface
# -------------------------------------------------------------------------
png(file.path(out_dir, "background_points_check.png"),
    width = 800, height = 800, res = 120)

par(mar = c(3, 3, 3, 5))
plot(bias_masked, main = "Background points over sampling bias surface")
points(bg_points$lon, bg_points$lat, pch = 20, cex = 0.3, col = "black")

dev.off()

# -------------------------------------------------------------------------
# 6. Save background points for ENMeval model-fitting
# -------------------------------------------------------------------------
out_bg_path <- file.path(base_dir, "background_points.csv")
write.csv(bg_points, out_bg_path, row.names = FALSE)

cat("Background points sampled, checked, and saved to:", out_bg_path, "\n")

# -------------------------------------------------------------------------
# 7. Quantitative check: did weighting actually bias the sample?
# -------------------------------------------------------------------------
sampled_weights <- bias_df$weight[bg_idx]
all_weights     <- bias_df$weight

cat("\n=== Diagnostic Weighting Metrics ===\n")
cat("Mean weight — all candidate cells:", mean(all_weights), "\n")
cat("Mean weight — sampled background points:", mean(sampled_weights), "\n")
cat("Ratio (should be > 1 if weighting worked):",
    round(mean(sampled_weights) / mean(all_weights), 3), "\n")

# Save side-by-side histograms
hist_path <- file.path(base_dir, "bias_weight_histograms.png")
png(hist_path, width = 1000, height = 500, res = 120)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(all_weights, main = "All Candidate Cells", xlab = "Bias Weight", col = "skyblue", breaks = 30)
hist(sampled_weights, main = "Sampled Background Points", xlab = "Bias Weight", col = "lightgreen", breaks = 30)

dev.off()

# Reset graphics parameters
par(mfrow = c(1, 1))

cat("Histograms successfully saved as 'bias_weight_histograms.png'!\n")