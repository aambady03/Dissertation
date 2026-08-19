# =========================================================
# Step 2: Generate Sampling Bias Surface for MaxEnt
# =========================================================
library(terra)

# -------------------------------------------------------------------------
# Load Environmental Layers & Calibration Area
# -------------------------------------------------------------------------

climate_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/climate/processed/"
climate_files <- list.files(climate_dir, pattern = "\\.tif$", full.names = TRUE)
env           <- rast(climate_files)

shapefile_path <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/calibration_area.shp"
m_area         <- vect(shapefile_path)

# -------------------------------------------------------------------------
# Quick Variable Check
# -------------------------------------------------------------------------
print("--- Environmental Layer Names ---")
print(names(env))

print("--- Coordinate Reference System Alignment ---")
print(crs(env, describe = TRUE)$name)
print(crs(m_area, describe = TRUE)$name)

# -------------------------------------------------------------------------
# 1. Load Local Major Roads Shapefile
# -------------------------------------------------------------------------
print("Loading local major roads dataset...")

roads_path <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/ne_10m_roads/ne_10m_roads.shp"

if (!file.exists(roads_path)) {
  stop("Please download the Natural Earth roads shapefile and save it to: ", roads_path)
}

global_roads <- vect(roads_path)

print("Cropping roads to calibration area...")
roads <- crop(global_roads, m_area)
roads <- project(roads, crs(env))

# -------------------------------------------------------------------------
# 2. Rasterize to climate grid
# -------------------------------------------------------------------------
road_rast <- rasterize(roads, env[["mean_ann_temp"]], field = 1, background = 0)

# -------------------------------------------------------------------------
# 3. Kernel Density Smoothing (Gaussian Focal Filter)
# -------------------------------------------------------------------------
print("Applying spatial Gaussian smoothing filter...")

# d = 0.5 degrees (~50km sigma). focalMat creates a 2D Gaussian kernel.
fw <- focalMat(road_rast, d = 0.5, type = "Gauss")
bias_file <- focal(road_rast, w = fw, fun = sum, na.rm = TRUE)

# Clip cleanly to your land boundary
bias_file <- crop(bias_file, m_area) |> mask(m_area)

# -------------------------------------------------------------------------
# 4. Correct Normalization to strict [0, 1] range
# -------------------------------------------------------------------------
# Shift values so minimum is non-zero (1e-4), then scale maximum to 1.0
b_min <- minmax(bias_file)[1]

# Add tiny constant to non-sampled areas to avoid absolute zero barriers
bias_file_offset <- (bias_file - b_min) + 0.0001
b_max <- minmax(bias_file_offset)[2]

bias_file_norm <- bias_file_offset / b_max

# -------------------------------------------------------------------------
# 5. Export for MaxEnt
# -------------------------------------------------------------------------
output_file <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/bias/results/bias_file.tif"
output_dir  <- dirname(output_file)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

writeRaster(bias_file_norm, output_file, overwrite = TRUE)
print("Successfully generated and saved your sampling bias surface!")

# Diagnostic Plot
png("/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/bias_surface_plot.png", 
    width = 800, height = 800, res = 120)
par(mar = c(3, 3, 3, 5))
plot(bias_file_norm, main = "Sampling Bias Surface (0.0001 - 1.0)")
dev.off()

print("Plot successfully saved directly to your folder as a PNG file!")