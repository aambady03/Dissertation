# =========================================================
# Step 1: Occurrence Thinning & Calibration Zone (M)
# =========================================================

library(spThin)
library(terra)
library(geodata)

# Define directories
base_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"
thin_dir <- file.path(base_dir, "thinned_data")

if (!dir.exists(thin_dir)) dir.create(thin_dir, recursive = TRUE)

# ---------------------------------------------------------
# 1. Load GBIF Data & Spatial Thinning
# ---------------------------------------------------------
gbif <- read.csv(file.path(base_dir, "parthenium_india_gbif.csv"))

# Quick clean: remove missing coordinates
clean <- gbif[!is.na(gbif$decimalLatitude) & !is.na(gbif$decimalLongitude), ]

# Spatial thinning (5 km distance threshold)
thinned <- spThin::thin(
  loc.data    = clean,
  lat.col     = "decimalLatitude",
  long.col    = "decimalLongitude",
  spec.col    = "species",
  thin.par    = 5,
  reps        = 10,
  out.dir     = thin_dir,
  out.base    = "parthenium_thinned",
  write.files = TRUE
)

# ---------------------------------------------------------
# Select & Load Best Thinning Replicate (Max Retained Points)
# ---------------------------------------------------------
# 1. Find all generated thinning CSV files in thin_dir
thin_files <- list.files(
  path       = thin_dir, 
  pattern    = "^parthenium_thinned_thin.*\\.csv$", 
  full.names = TRUE
)

# 2. Identify the file containing the maximum number of records
file_row_counts <- sapply(thin_files, function(f) nrow(read.csv(f)))
best_thin_file  <- thin_files[which.max(file_row_counts)]

cat(sprintf("Selected best replicate: %s (%d records retained)\n", 
            basename(best_thin_file), max(file_row_counts)))

# 3. Load the best replicate
occs_raw <- read.csv(best_thin_file)
occs_raw <- occs_raw[, c("decimalLongitude", "decimalLatitude")]
names(occs_raw) <- c("lon", "lat")

cat("Thinned records before land boundary filtering:", nrow(occs_raw), "\n")

# ---------------------------------------------------------
# 2. Land Boundary Filtering & Calibration Area (M)
# ---------------------------------------------------------

# Convert raw occurrence points to SpatVector
pts_raw <- vect(occs_raw, geom = c("lon", "lat"), crs = "EPSG:4326")

# Load India administrative boundary map
india <- gadm("IND", level = 0, path = base_dir)

# Filter out occurrence points falling outside land borders
pts_clean <- pts_raw[india, ]

# Extract clean spatial points dataframe
occs_clean <- as.data.frame(pts_clean, geom = "XY")
colnames(occs_clean) <- c("lon", "lat")
cat("Cleaned records within land boundaries:", nrow(occs_clean), "\n")

# CRUCIAL FIX: Write cleaned occurrences to disk for downstream scripts
clean_occ_path <- file.path(base_dir, "parthenium_thinned_clean.csv")
write.csv(occs_clean, clean_occ_path, row.names = FALSE)
cat("Saved cleaned occurrences to:", clean_occ_path, "\n")

# ---------------------------------------------------------
# 3. Buffer Construction (100km M Area)
# ---------------------------------------------------------

# Project points to metric system (UTM Zone 43N / EPSG:7755) for accurate buffering
pts_metric <- project(pts_clean, "EPSG:7755")

# Create 100 km buffers around clean occurrence points
buffered_pts_metric <- buffer(pts_metric, width = 100000)

# Re-project buffers back to WGS84 (EPSG:4326)
m_area_raw <- project(buffered_pts_metric, "EPSG:4326")

# Dissolve overlapping circles and fill interior gaps/holes
m_area_dissolved <- aggregate(m_area_raw)
m_area_filled    <- fillHoles(m_area_dissolved)

# Clip buffer boundary cleanly to India land borders
calibration_area <- intersect(m_area_filled, india)

# Save calibration area shapefile
calib_out_path <- file.path(base_dir, "calibration_area.shp")
writeVector(calibration_area, calib_out_path, overwrite = TRUE)
cat("Calibration area (M) saved to:", calib_out_path, "\n")

# ---------------------------------------------------------
# 4. Diagnostic Visual Check
# ---------------------------------------------------------
png(file.path(base_dir, "calibration_area_check.png"), 
    width = 800, height = 800, res = 120)

plot(india, main = "Cleaned Calibration Area (M) & Presences")
plot(calibration_area, add = TRUE, border = "blue", col = rgb(0, 0, 1, 0.15))
points(pts_clean, cex = 0.4, col = "red", pch = 16)

dev.off()

cat("\nStep 1 execution complete! 'parthenium_thinned_clean.csv' and 'calibration_area.shp' ready.\n")

# =========================================================
# 5. Summary counts for Methods write-up
# =========================================================

n_gbif_total <- nrow(gbif)
n_with_coords <- nrow(clean)
n_after_thin <- nrow(occs_raw)
n_after_land_filter <- nrow(occs_clean)
n_thin_reps <- 10
thin_distance_km <- 5

cat("\n--- Occurrence data summary ---\n")
cat("GBIF records downloaded:", n_gbif_total, "\n")
cat("Records with valid coordinates:", n_with_coords, "\n")
cat("Spatial thinning distance (km):", thin_distance_km, "\n")
cat("Number of thinning replicates:", n_thin_reps, "\n")
cat("Records retained in selected thinning replicate:", n_after_thin, "\n")
cat("Records retained after land-boundary filtering:", n_after_land_filter, "\n")
cat("Selected thinning replicate:", basename(best_thin_file), "\n")

summary_df <- data.frame(
  stage = c(
    "GBIF downloaded",
    "Valid coordinates retained",
    "Selected thinning replicate retained",
    "After India land-boundary filtering"
  ),
  n = c(
    nrow(gbif),
    nrow(clean),
    nrow(occs_raw),
    nrow(occs_clean)
  )
)

write.csv(
  summary_df,
  file.path(base_dir, "occurrence_data_summary.csv"),
  row.names = FALSE
)

cat("Saved occurrence summary to:", file.path(base_dir, "occurrence_data_summary.csv"), "\n")