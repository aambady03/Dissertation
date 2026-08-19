# =========================================================
# Step 5: MaxEnt Model Fitting & Evaluation (ENMeval)
# =========================================================

library(ENMeval)
library(terra)
library(maxnet)

base_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"

# -------------------------------------------------------------------------
# 1. Load Environmental Predictors
# -------------------------------------------------------------------------
climate_dir <- file.path(base_dir, "climate/processed")

# Explicitly load rasters to guarantee layer naming and ordering
mean_ann_temp    <- rast(file.path(climate_dir, "mean_ann_temp.tif"))
total_ann_precip <- rast(file.path(climate_dir, "total_ann_precip.tif"))
temp_seasonality <- rast(file.path(climate_dir, "temp_seasonality.tif"))

envs <- c(mean_ann_temp, total_ann_precip, temp_seasonality)
names(envs) <- c("mean_ann_temp", "total_ann_precip", "temp_seasonality")

cat("Layer names:", names(envs), "\n")
cat("Layer order confirmed — mean_ann_temp | total_ann_precip | temp_seasonality\n")
print(envs)

# -------------------------------------------------------------------------
# 2. Load Occurrence and Background Points
# -------------------------------------------------------------------------
occs_path <- file.path(base_dir, "parthenium_thinned_clean.csv")
bg_path   <- file.path(base_dir, "background_points.csv")

if (!file.exists(occs_path) || !file.exists(bg_path)) {
  stop("Missing thinned occurrences or sampled background points. Run Steps 1 & 3 first!")
}

occs <- read.csv(occs_path)
# Handle variations in longitude/latitude column naming
if (all(c("decimalLongitude", "decimalLatitude") %in% names(occs))) {
  occs <- occs[, c("decimalLongitude", "decimalLatitude")]
} else {
  occs <- occs[, c("lon", "lat")]
}
names(occs) <- c("lon", "lat")

bg_points <- read.csv(bg_path)
if (all(c("decimalLongitude", "decimalLatitude") %in% names(bg_points))) {
  bg_points <- bg_points[, c("decimalLongitude", "decimalLatitude")]
  names(bg_points) <- c("lon", "lat")
}

cat("Raw occurrence points loaded:", nrow(occs), "\n")
cat("Background points loaded:", nrow(bg_points), "\n")

# -------------------------------------------------------------------------
# 3. Filter NA Extraction Cells & Save Clean Sets
# -------------------------------------------------------------------------
occs_vals <- terra::extract(envs, occs[, c("lon", "lat")])
bg_vals   <- terra::extract(envs, bg_points[, c("lon", "lat")])

occs_na <- which(!complete.cases(occs_vals))
bg_na   <- which(!complete.cases(bg_vals))

if (length(occs_na) > 0) {
  cat("Dropping", length(occs_na), "occurrence points falling outside raster grid / NA cells\n")
  occs <- occs[-occs_na, ]
}
if (length(bg_na) > 0) {
  cat("Dropping", length(bg_na), "background points falling outside raster grid / NA cells\n")
  bg_points <- bg_points[-bg_na, ]
}

cat("Final valid occurrence points:", nrow(occs), "\n")
cat("Final valid background points:", nrow(bg_points), "\n")

# Save cleaned occurrences for downstream MESS reference
write.csv(occs, file.path(base_dir, "occs_clean.csv"), row.names = FALSE)
cat("Cleaned occurrences written to occs_clean.csv\n")

# -------------------------------------------------------------------------
# 4. Tune MaxEnt via ENMeval
# -------------------------------------------------------------------------
e <- ENMevaluate(
  occs       = occs,
  envs       = envs,
  bg         = bg_points,
  algorithm  = "maxnet",
  partitions = "block",
  tune.args  = list(fc = c("L", "LQ", "LQH"), rm = seq(0.5, 4, 0.5))
)

# -------------------------------------------------------------------------
# 5. Select Optimal Model using AICc (Delta AICc = 0)
# -------------------------------------------------------------------------
results_tbl <- eval.results(e)

# Sort by delta.AICc ascending and select top model
best_row <- results_tbl[order(results_tbl$delta.AICc), ][1, ]

cat("\n=== BEST MODEL SELECTION (AICc) ===\n")
cat("Selected Model Code :", best_row$tune.args, "\n")
cat("AICc Value          :", round(best_row$AICc, 3), "\n")
cat("Delta AICc          :", round(best_row$delta.AICc, 3), "\n")
cat("Parameters (K)      :", best_row$ncoef, "\n\n")

mod_name <- best_row$tune.args
best_mod <- eval.models(e)[[mod_name]]

# -------------------------------------------------------------------------
# 6. Predict Current Suitability & Calculate Thresholds
# -------------------------------------------------------------------------
pred_india <- predict(envs, best_mod, type = "cloglog", na.rm = TRUE)

raster_min <- minmax(pred_india)[1]
raster_max <- minmax(pred_india)[2]
cat(sprintf("Suitability Range: %.3f to %.3f\n", raster_min, raster_max))

if (raster_max > 1 || raster_min < 0) {
  stop("Suitability values fall outside the valid 0–1 cloglog bounds!")
}

# Extract values at presence locations
pred_vals <- terra::extract(pred_india, occs[, c("lon", "lat")])
suitability_at_occs <- pred_vals[[2]]

cat("Mean suitability at presence points  :", round(mean(suitability_at_occs, na.rm = TRUE), 3), "\n")
cat("Median suitability at presence points:", round(median(suitability_at_occs, na.rm = TRUE), 3), "\n")

# Calculate 10th percentile training presence (P10) threshold
p10_threshold <- quantile(suitability_at_occs, 0.10, na.rm = TRUE)
cat("P10 Threshold (10th Percentile)      :", round(p10_threshold, 4), "\n")

saveRDS(p10_threshold, file.path(base_dir, "p10_threshold.rds"))

# Binarize current suitability map
suit_binary <- pred_india >= p10_threshold

# Save Raster Outputs
writeRaster(pred_india, file.path(base_dir, "parthenium_suitability.tif"), overwrite = TRUE)
writeRaster(suit_binary, file.path(base_dir, "parthenium_suitability_binary.tif"), overwrite = TRUE)

# Diagnostic Plots
png(file.path(base_dir, "parthenium_suitability_maps.png"), width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
plot(pred_india, main = "Parthenium Suitability (Cloglog)")
plot(suit_binary, main = "Binary Suitability (>= P10 Threshold)")
par(mfrow = c(1, 1))
dev.off()

# -------------------------------------------------------------------------
# 7. Save Tuning Results & Model Object
# -------------------------------------------------------------------------
write.csv(results_tbl, file.path(base_dir, "enmeval_results_table.csv"), row.names = FALSE)
saveRDS(best_mod, file.path(base_dir, "best_maxnet_model.rds"))

cat("\nStep 5 Model Training Complete! Best model object saved as 'best_maxnet_model.rds'.\n")