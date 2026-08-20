library(geodata)
library(terra)

# 1. Load baseline reference layer
mat_curr <- rast("/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/climate/processed/mean_ann_temp.tif")

# 2. Download CMIP6 world layers (CanESM5, SSP2-4.5)
cat("Downloading raw CMIP6 bioclimatic data...\n")
bio_2030_raw <- cmip6_world(var = "bio", model = "CanESM5", ssp = "245", time = "2021-2040", res = 2.5, path = tempdir())
bio_2050_raw <- cmip6_world(var = "bio", model = "CanESM5", ssp = "245", time = "2041-2060", res = 2.5, path = tempdir())

# 3. Extract BIO5 (Max Temperature of Warmest Month = Layer 5)
summer_heat_2030 <- bio_2030_raw[[5]]
summer_heat_2050 <- bio_2050_raw[[5]]

names(summer_heat_2030) <- "max_summer_temp"
names(summer_heat_2050) <- "max_summer_temp"

# 4. Crop & Mask strictly to India spatial extent using your reference layer
cat("Cropping and aligning to India spatial extent...\n")
summer_heat_2030 <- crop(summer_heat_2030, mat_curr) |> resample(mat_curr) |> mask(mat_curr)
summer_heat_2050 <- crop(summer_heat_2050, mat_curr) |> resample(mat_curr) |> mask(mat_curr)

# 5. Save extreme summer heat rasters
path_2030 <- "~/Documents/DISSERT/EDITS/section 2/climate/future/max_summer_temp_2030.tif"
path_2050 <- "~/Documents/DISSERT/EDITS/section 2/climate/future/max_summer_temp_2050.tif"

writeRaster(summer_heat_2030, path_2030, overwrite = TRUE)
writeRaster(summer_heat_2050, path_2050, overwrite = TRUE)

cat("Successfully saved:\n - ", path_2030, "\n - ", path_2050, "\n")