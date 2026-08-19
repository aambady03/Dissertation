library(terra)
library(ncdf4)

# =========================================================
# 1. Load Calibration Area (M) from Step 1
# =========================================================
base_dir   <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"
calib_path <- file.path(base_dir, "calibration_area.shp")

if (!file.exists(calib_path)) {
  stop("calibration_area.shp not found. Run calibration_zone.R first!")
}

calibration_area <- vect(calib_path)

# =========================================================
# 2. Process ERA5-Land Monthly Climatology (1991–2020)
# =========================================================
nc_path <- file.path(base_dir, "climate/era5_land_monthly_1991_2020.nc")

if (!file.exists(nc_path)) {
  stop("ERA5 NetCDF file not found at: ", nc_path)
}

# Load multi-variable raster stack from NetCDF
climate_raw <- rast(nc_path)

# Ensure CRS is explicitly WGS84 (EPSG:4326) to match shapefiles
if (is.na(crs(climate_raw)) || crs(climate_raw) == "") {
  crs(climate_raw) <- "EPSG:4326"
}

# Isolate temperature (t2m) and precipitation (tp) layers
temp_layers   <- climate_raw[[grep("t2m", names(climate_raw))]]
precip_layers <- climate_raw[[grep("tp",  names(climate_raw))]]

# =========================================================
# 3. Derive Long-Term Climatological Summaries
# =========================================================

# 1. Mean Annual Temperature: Convert Kelvin to Celsius
mat <- mean(temp_layers) - 273.15  

# 2. Total Annual Precipitation: 
# ERA5 tp is in m/day. 
# Mean daily rate over 30 yrs (m/day) * 365.25 days/yr * 1000 mm/m = mm/year
mean_daily_rate <- mean(precip_layers) # mean m/day
tap <- mean_daily_rate * 365.25 * 1000  

# 3. Temperature Seasonality: Standard deviation across monthly layers in °C
seas <- stdev(temp_layers)         

# Combine into single environmental stack
env_full <- c(mat, tap, seas)
names(env_full) <- c("mean_ann_temp", "total_ann_precip", "temp_seasonality")

# =========================================================
# 4. Spatial Masking & Export
# =========================================================

# Ensure CRS alignment between raster stack and vector polygon
if (crs(env_full) != crs(calibration_area)) {
  calibration_area <- project(calibration_area, crs(env_full))
}

# Crop and mask directly to calibration area
env <- crop(env_full, calibration_area) |> mask(calibration_area)

# Create output directory if it doesn't exist
out_dir <- file.path(base_dir, "climate/processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save individual layer rasters for downstream scripts
writeRaster(env[["mean_ann_temp"]],    file.path(out_dir, "mean_ann_temp.tif"),    overwrite = TRUE)
writeRaster(env[["total_ann_precip"]],  file.path(out_dir, "total_ann_precip.tif"),  overwrite = TRUE)
writeRaster(env[["temp_seasonality"]], file.path(out_dir, "temp_seasonality.tif"), overwrite = TRUE)

# Diagnostic verification
cat("\n=== Sanity Check on Processed Climate Ranges ===\n")
print(summary(env))

# Save Diagnostic plot
png(file.path(base_dir, "climate_layers_check.png"), 
    width = 1000, height = 400, res = 120)
plot(env)
dev.off()

cat("\nClimate variables successfully processed, masked, and saved!\n")