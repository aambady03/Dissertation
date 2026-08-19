# ==============================================================================
# SCRIPT 01: DATA PREPARATION & INGESTION
# ==============================================================================

library(tidyverse)

# 1. Load Raw CSV (edit to data folder)
raw_df <- read.csv("/home/anaga-ambady/Documents/DISSERT/DATA/python/results/parthenium_with_climate_matched.csv")

# 2. Process and Clean Data
survey_data <- raw_df %>%
  mutate(
    hab_group = as.factor(case_when(
      hab %in% c("Field crop - planted", "Orchard")                              ~ "Agricultural_managed",
      hab %in% c("Field crop- margin", "Field crop - fallow", "Field crop - bare soil") ~ "Agricultural_disturbed",
      hab %in% c("Grassland", "Shrubland and woodland", "Bare areas")           ~ "Semi_natural_open",
      hab %in% c("Building, roads or runways", "Parks, gardens and cemeteries")  ~ "Urban_disturbed",
      hab %in% c("Lakes, ponds, rivers and streams", "Natural vegetation in wet areas") ~ "Wet_habitat",
      hab %in% c("Forest", "Agroforestry")                                     ~ "Woody",
      TRUE ~ "Other"
    )),
    road_type  = as.factor(road_type),
    plotID     = as.factor(plotID),
    longitude  = as.numeric(longitude),
    latitude   = as.numeric(latitude)
  ) %>%
  filter(
    !is.na(presence),
    !is.na(longitude),
    !is.na(latitude),
    !is.na(plotID),
    hab_group != "Other"
  )

# 3. Validation Prints
cat("Total Observations:", nrow(survey_data), "\n")
cat("Presences:", sum(survey_data$presence == 1), "\n")
cat("Absences:", sum(survey_data$presence == 0), "\n")
print(table(survey_data$hab_group))

# 4. Save Clean Dataset
out_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/MODEL/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(survey_data, file.path(out_dir, "survey_data.csv"), row.names = FALSE)