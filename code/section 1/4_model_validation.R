# ==============================================================================
# SCRIPT 04: COMPREHENSIVE DIAGNOSTICS & VALIDATION
# ==============================================================================

library(tidyverse)
library(mgcv)
library(DHARMa)
library(pROC)
library(caret)
library(rms)
library(ecospat)
library(sf)
library(spdep)
library(blockCV)

out_dir     <- "~/Documents/DISSERT/EDITS/section 1/results"
survey_data <- read.csv("~/Documents/DISSERT/EDITS/section 1/survey_data.csv")
gam_gdd_k_200 <- readRDS("~/Documents/DISSERT/EDITS/section 1/gam_rds/gam_gdd_k_200.rds")

# --- 0. Ensure Projected Coordinates (x_km, y_km) Exist ---
if (!all(c("x_km", "y_km") %in% colnames(survey_data))) {
  pts_proj <- st_as_sf(survey_data, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(32643) # UTM Zone 43N
  coords_m <- st_coordinates(pts_proj)
  survey_data$x_km <- coords_m[, 1] / 1000
  survey_data$y_km <- coords_m[, 2] / 1000
}

# --- 1. Concurvity Diagnostics ---
conc_full <- concurvity(gam_gdd_k_200, full = TRUE)
conc_pair <- concurvity(gam_gdd_k_200, full = FALSE)
write.csv(conc_full, file.path(out_dir, "concurvity_full.csv"))
write.csv(conc_pair, file.path(out_dir, "concurvity_pairwise.csv"))

# --- 2. Matched-k Altitude Test ---
gam_no_alt_matched <- update(
  gam_gdd_k_200, 
  . ~ . - s(altitude), 
  data = gam_gdd_k_200$model
)

# --- 3. Spatial Autocorrelation (Moran's I & DHARMa) ---
survey_data$resid <- residuals(gam_gdd_k_200, type = "pearson")

# Raw Pearson Moran's I using metric coordinates (meters)
coords_m <- cbind(survey_data$x_km * 1000, survey_data$y_km * 1000)
coords_m_jittered <- st_coordinates(st_jitter(st_as_sf(as.data.frame(coords_m), coords = c(1, 2)), amount = 1))

moran_raw <- moran.test(survey_data$resid, nb2listw(knn2nb(knearneigh(coords_m_jittered, k = 8))))

write.csv(
  data.frame(Moran_I = moran_raw$estimate["Moran I statistic"], p_val = moran_raw$p.value), 
  file.path(out_dir, "moran_raw_pearson.csv"), 
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Aggregated DHARMa Moran's I Test & Plotting (Using x_km, y_km)
# ------------------------------------------------------------------------------
set.seed(42)
sim_res <- simulateResiduals(gam_gdd_k_200, n = 500)

# Group residuals by spatial location (x_km, y_km)
sim_res_agg <- recalculateResiduals(
  sim_res, 
  group = interaction(survey_data$x_km, survey_data$y_km, drop = TRUE)
)

# Extract matching unique coordinates in exact grouped order
unique_coords <- survey_data %>%
  dplyr::group_by(x_km, y_km) %>%
  dplyr::summarise(.groups = "drop")

# 1. Run DHARMa Spatial Autocorrelation Test
moran_dharma <- testSpatialAutocorrelation(
  sim_res_agg, 
  x = unique_coords$x_km, 
  y = unique_coords$y_km,
  plot = FALSE
)

write.csv(
  data.frame(Moran_I = moran_dharma$statistic, p_val = moran_dharma$p.value), 
  file.path(out_dir, "moran_dharma_aggregated.csv"), 
  row.names = FALSE
)

# 2. Save DHARMa spatial plot to disk safely
png(file.path(out_dir, "dharma_spatial_autocorrelation.png"), width = 800, height = 600)
plotResiduals(sim_res_agg, form = unique_coords$x_km)
dev.off()

# --- 4. Extreme Residual Analysis ---
extreme_rows <- which(survey_data$resid > 10)
for(row_i in extreme_rows) {
  nb <- survey_data[abs(survey_data$x_km - survey_data$x_km[row_i]) < 1 & 
                      abs(survey_data$y_km - survey_data$y_km[row_i]) < 1, ]
  write.csv(nb, file.path(out_dir, paste0("neighbourhood_row_", row_i, ".csv")), row.names = FALSE)
}

# --- 5. ROC / AUC Discrimination & Confusion Matrix ---
pred_prob <- predict(gam_gdd_k_200, type = "response")
roc_obj   <- roc(survey_data$presence, pred_prob)
opt_thresh <- roc_obj$thresholds[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)]

cm_youden <- confusionMatrix(factor(ifelse(pred_prob >= opt_thresh, 1, 0), levels = c(0,1)), factor(survey_data$presence, levels = c(0,1)), positive = "1")
cm_05     <- confusionMatrix(factor(ifelse(pred_prob >= 0.5, 1, 0), levels = c(0,1)), factor(survey_data$presence, levels = c(0,1)), positive = "1")

thresh_df <- data.frame(
  Threshold = c("Youden J", "0.5 Fixed"), Value = c(opt_thresh, 0.5),
  Sensitivity = c(cm_youden$byClass["Sensitivity"], cm_05$byClass["Sensitivity"]),
  Specificity = c(cm_youden$byClass["Specificity"], cm_05$byClass["Specificity"]),
  Kappa = c(cm_youden$overall["Kappa"], cm_05$overall["Kappa"])
)
write.csv(thresh_df, file.path(out_dir, "threshold_comparison.csv"), row.names = FALSE)

# --- 6. Calibration & Boyce Index ---
png(file.path(out_dir, "calibration_plot.png"), width = 800, height = 600)
val.prob(pred_prob, survey_data$presence, pl = TRUE, smooth = TRUE)
dev.off()

png(file.path(out_dir, "boyce_index_plot.png"), width = 800, height = 600)
boyce_full <- ecospat.boyce(fit = pred_prob, obs = pred_prob[survey_data$presence == 1], nclass = 10, PEplot = TRUE)
dev.off()

# ==============================================================================
# 7. Spatial Block Cross-Validation (50 km) - ROBUST FOLD HANDLER
# ==============================================================================

survey_sf <- st_as_sf(survey_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Generate 50 km spatial blocks
cv_blocks <- cv_spatial(x = survey_sf, column = "presence", k = 5, size = 50000)

auc_cv <- sens_cv <- spec_cv <- boyce_cv <- numeric(length(cv_blocks$folds_list))

for (i in seq_along(cv_blocks$folds_list)) {
  cat(sprintf("\n--- Processing Fold %d of %d ---\n", i, length(cv_blocks$folds_list)))
  
  # 1. Extract train & test sets and drop SF geometry
  train_df <- st_drop_geometry(survey_sf[cv_blocks$folds_list[[i]][[1]], ])
  test_df  <- st_drop_geometry(survey_sf[cv_blocks$folds_list[[i]][[2]], ])
  
  # 2. Filter out any missing values in predictors for this fold
  model_vars <- c("presence", "hab_group", "road_type", "precip_sum_30d", 
                  "gdd_30d", "altitude", "hab_perc", "x_km", "y_km")
  
  train_df <- train_df %>% filter(if_all(all_of(model_vars), ~ !is.na(.)))
  test_df  <- test_df  %>% filter(if_all(all_of(model_vars), ~ !is.na(.)))
  
  # 3. Handle factor levels cleanly
  train_df$hab_group <- factor(train_df$hab_group, levels = levels(survey_data$hab_group))
  train_df$road_type <- factor(train_df$road_type, levels = levels(survey_data$road_type))
  test_df$hab_group  <- factor(test_df$hab_group,  levels = levels(survey_data$hab_group))
  test_df$road_type  <- factor(test_df$road_type,  levels = levels(survey_data$road_type))
  
  # 4. Dynamic Basis Dimension (k) Check for Small Folds
  n_train <- nrow(train_df)
  cat(sprintf("Training observations in Fold %d: %d\n", i, n_train))
  
  if (n_train < 100) {
    warning(sprintf("Fold %d has only %d rows. Reducing smooth basis dimensions (k).", i, n_train))
  }
  
  # Cap k dynamically based on fold size to avoid "Not enough data" error
  k_sp    <- min(80, max(10, floor(n_train * 0.15)))
  k_clim  <- min(25, max(5,  floor(n_train * 0.05)))
  k_hab   <- min(15, max(4,  floor(n_train * 0.03)))
  
  # 5. Fit Fold Model Matching M5 Specification
  fmod <- tryCatch({
    gam(
      presence ~ hab_group + road_type + 
        s(precip_sum_30d, k = k_clim) + 
        s(gdd_30d, k = k_clim) +
        s(gdd_30d, by = hab_group, k = k_hab) + 
        s(altitude, k = k_hab) + 
        s(hab_perc, k = 10) +
        s(x_km, y_km, bs = "gp", k = k_sp),
      family = binomial(link = "logit"), 
      method = "REML", 
      data   = train_df
    )
  }, error = function(e) {
    message(sprintf("Error fitting fold %d: %s", i, e$message))
    return(NULL)
  })
  
  # Skip evaluation if fold model failed
  if (is.null(fmod)) {
    auc_cv[i]   <- NA
    sens_cv[i]  <- NA
    spec_cv[i]  <- NA
    boyce_cv[i] <- NA
    next
  }
  
  # 6. Predict on test fold
  p_fold <- predict(fmod, newdata = test_df, type = "response")
  
  # Metrics
  auc_cv[i] <- as.numeric(auc(test_df$presence, p_fold, quiet = TRUE))
  
  tbl <- table(
    factor(ifelse(p_fold >= 0.5, 1, 0), levels = c(0,1)), 
    factor(test_df$presence, levels = c(0,1))
  )
  sens_cv[i] <- tbl[2,2] / sum(tbl[,2])
  spec_cv[i] <- tbl[1,1] / sum(tbl[,1])
  
  boyce_cv[i] <- tryCatch({
    ecospat.boyce(fit = p_fold, obs = p_fold[test_df$presence == 1], nclass = 10, PEplot = FALSE)$Spearman.cor
  }, error = function(e) NA)
}

# Summary output
cv_summary <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "Boyce Index"),
  Mean   = round(c(mean(auc_cv, na.rm=TRUE), mean(sens_cv, na.rm=TRUE), mean(spec_cv, na.rm=TRUE), mean(boyce_cv, na.rm = TRUE)), 3),
  SD     = round(c(sd(auc_cv, na.rm=TRUE), sd(sens_cv, na.rm=TRUE), sd(spec_cv, na.rm=TRUE), sd(boyce_cv, na.rm = TRUE)), 3)
)

print(cv_summary)
write.csv(cv_summary, file.path(out_dir, "spatial_block_cv_summary.csv"), row.names = FALSE)