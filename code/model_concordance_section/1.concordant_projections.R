## =========================================================================
## IMPROVED WORKFLOW: GAM & MaxEnt Spatial Projections + Concordance Analysis
## =========================================================================
## Purpose:
##   1. Project India-calibrated MaxEnt across Pakistan
##   2. Project Pakistan-calibrated GAM across Pakistan
##   3. Force identical spatial support and identical valid-cell comparison set
##   4. Compute correlation, overlap, concordance, and rank-difference metrics
##   5. Produce a publication-ready 4-panel figure with colour-blind friendly palettes
## =========================================================================

# -------------------------------------------------------------------------
# STEP 0: Load packages
# -------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(mgcv)
  library(terra)
  library(maxnet)
  library(dplyr)
  library(ggplot2)
  library(tidyterra)
  library(patchwork)
  library(grid)
})

# -------------------------------------------------------------------------
# STEP 1: Paths and output directory
# -------------------------------------------------------------------------
base_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/bridge"
out_dir  <- "~/Documents/DISSERT/EDITS/bridge/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

maxent_model_path <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india/best_maxnet_model.rds"
gam_model_path    <- "~/Documents/DISSERT/EDITS/section 1/gam_rds/gam_gdd_k_200.rds"

# -------------------------------------------------------------------------
# STEP 2: Predict MaxEnt suitability across Pakistan
# -------------------------------------------------------------------------
message("--> Step 2: Running MaxEnt projection...")
best_mod <- readRDS(maxent_model_path)

pak_envs <- rast(c(
  file.path(base_dir, "pakistan_mean_ann_temp.tif"),
  file.path(base_dir, "pakistan_total_ann_precip.tif"),
  file.path(base_dir, "pakistan_temp_seasonality.tif")
))
names(pak_envs) <- c("mean_ann_temp", "total_ann_precip", "temp_seasonality")

maxent_pred_raw_r <- predict(pak_envs, best_mod, type = "cloglog", na.rm = TRUE)
writeRaster(
  maxent_pred_raw_r,
  file.path(out_dir, "maxent_pakistan_cloglog_raw.tif"),
  overwrite = TRUE
)
message("✓ Raw MaxEnt projection saved.")

# -------------------------------------------------------------------------
# STEP 3: Load GAM model and Pakistan rasters
# -------------------------------------------------------------------------
message("--> Step 3: Preparing GAM covariates...")
gam_model <- readRDS(gam_model_path)

gdd_r     <- rast(file.path(base_dir, "rasters", "pakistan_gdd_30d.tif"))
precip_r  <- rast(file.path(base_dir, "rasters", "pakistan_precip_30d.tif"))
alt_r     <- rast(file.path(base_dir, "rasters", "pakistan_altitude.tif"))
hab_pct_r <- rast(file.path(base_dir, "rasters", "pakistan_habitat_pct_cover.tif"))

# Match all continuous GAM covariates to the GDD reference grid
precip_r  <- resample(precip_r,  gdd_r, method = "bilinear")
alt_r     <- resample(alt_r,     gdd_r, method = "bilinear")
hab_pct_r <- resample(hab_pct_r, gdd_r, method = "bilinear")

# -------------------------------------------------------------------------
# STEP 4: Assign fixed factor reference levels from GAM training data
# -------------------------------------------------------------------------
hab_group_ref <- names(sort(table(gam_model$model$hab_group), decreasing = TRUE))[1]
road_type_ref <- names(sort(table(gam_model$model$road_type), decreasing = TRUE))[1]

message("hab_group fixed to modal level: ", hab_group_ref)
message("road_type fixed to modal level: ", road_type_ref)

# -------------------------------------------------------------------------
# STEP 5: Build prediction dataframe and enforce environmental envelope
# -------------------------------------------------------------------------
covariate_stack <- c(gdd_r, precip_r, alt_r, hab_pct_r)
names(covariate_stack) <- c("gdd_30d", "precip_sum_30d", "altitude", "hab_perc")

pred_df <- as.data.frame(covariate_stack, xy = TRUE, na.rm = TRUE)
names(pred_df)[names(pred_df) == "x"] <- "longitude"
names(pred_df)[names(pred_df) == "y"] <- "latitude"

# Add dummy spatial coordinates required by predict.gam (excluded during prediction)
pred_df$x_km <- 0
pred_df$y_km <- 0

# Assign fixed reference levels for factors
pred_df$hab_group <- factor(hab_group_ref, levels = levels(gam_model$model$hab_group))
pred_df$road_type <- factor(road_type_ref, levels = levels(gam_model$model$road_type))

# Enforce environmental envelope from model training range
mf <- gam_model$model
gdd_range    <- range(mf$gdd_30d, na.rm = TRUE)
precip_range <- range(mf$precip_sum_30d, na.rm = TRUE)
alt_range    <- range(mf$altitude, na.rm = TRUE)

pred_df <- pred_df %>%
  filter(
    gdd_30d        >= gdd_range[1]    & gdd_30d        <= gdd_range[2],
    precip_sum_30d >= precip_range[1] & precip_sum_30d <= precip_range[2],
    altitude       >= alt_range[1]    & altitude       <= alt_range[2]
  )

message("✓ Pixels inside training environmental envelope: ", format(nrow(pred_df), big.mark = ","))

# -------------------------------------------------------------------------
# STEP 6: Predict GAM and rasterise onto Pakistan grid
# -------------------------------------------------------------------------
message("--> Step 6: Predicting GAM suitability...")

# Exclude 2D spatial Gaussian Process smoother s(x_km, y_km) during regional transfer
pred_df$gam_suitability <- predict(
  gam_model,
  newdata = pred_df,
  type = "response",
  exclude = "s(x_km,y_km)"
)

gam_pred_r <- rast(
  pred_df[, c("longitude", "latitude", "gam_suitability")],
  type = "xyz",
  crs = crs(gdd_r)
)

writeRaster(
  gam_pred_r,
  file.path(out_dir, "gam_bridge_prediction.tif"),
  overwrite = TRUE
)
message("✓ GAM prediction raster saved.")

# -------------------------------------------------------------------------
# STEP 7: Align MaxEnt to GAM footprint explicitly
# -------------------------------------------------------------------------
message("--> Step 7: Aligning MaxEnt to GAM grid and footprint...")
maxent_pred_r <- resample(maxent_pred_raw_r, gam_pred_r, method = "bilinear")
maxent_pred_r <- mask(maxent_pred_r, gam_pred_r)

writeRaster(
  maxent_pred_r,
  file.path(out_dir, "maxent_pakistan_cloglog_aligned.tif"),
  overwrite = TRUE
)
message("✓ Aligned MaxEnt raster saved.")

# -------------------------------------------------------------------------
# STEP 8: Build one master comparison dataframe from overlapping valid cells
# -------------------------------------------------------------------------
message("--> Step 8: Building master comparison dataframe...")
gam_vals    <- values(gam_pred_r, mat = TRUE)[, 1]
maxent_vals <- values(maxent_pred_r, mat = TRUE)[, 1]
valid_cells <- which(!is.na(gam_vals) & !is.na(maxent_vals))

comparison_df <- data.frame(
  cell_id = valid_cells,
  gam     = gam_vals[valid_cells],
  maxent  = maxent_vals[valid_cells]
) %>%
  mutate(
    gam_rank    = percent_rank(gam),
    maxent_rank = percent_rank(maxent),
    rank_diff   = gam_rank - maxent_rank
  )

message("✓ Overlapping valid cells: ", format(nrow(comparison_df), big.mark = ","))

# -------------------------------------------------------------------------
# STEP 9: Top-quartile categorical agreement from the same master dataframe
# -------------------------------------------------------------------------
message("--> Step 9: Calculating categorical agreement...")
gam_q75    <- quantile(comparison_df$gam, 0.75, na.rm = TRUE)
maxent_q75 <- quantile(comparison_df$maxent, 0.75, na.rm = TRUE)

comparison_df <- comparison_df %>%
  mutate(
    gam_high    = gam > gam_q75,
    maxent_high = maxent > maxent_q75,
    agreement_code = case_when(
      !gam_high & !maxent_high ~ 0L,
      gam_high & !maxent_high  ~ 1L,
      !gam_high & maxent_high  ~ 2L,
      gam_high & maxent_high   ~ 3L
    )
  )

agreement_summary <- comparison_df %>%
  count(agreement_code, name = "count") %>%
  mutate(
    label = c("Both low", "GAM only high", "MaxEnt only high", "Both high")[agreement_code + 1],
    pct   = round(100 * count / sum(count), 1)
  )

cat("\n=== Top-Quartile Agreement Summary ===\n")
print(agreement_summary)

# -------------------------------------------------------------------------
# STEP 10: Reconstruct analysis rasters from the master dataframe
# -------------------------------------------------------------------------
message("--> Step 10: Rebuilding derived rasters from master dataframe...")
grid_template <- gam_pred_r
values(grid_template) <- NA

# Relative suitability maps as percentile ranks for better cross-model comparability
gam_rank_r <- grid_template
values(gam_rank_r)[comparison_df$cell_id] <- comparison_df$gam_rank

maxent_rank_r <- grid_template
values(maxent_rank_r)[comparison_df$cell_id] <- comparison_df$maxent_rank

agreement_r <- grid_template
values(agreement_r)[comparison_df$cell_id] <- comparison_df$agreement_code
agreement_factor <- as.factor(agreement_r)
levels(agreement_factor) <- data.frame(
  id    = 0:3,
  label = c("Both low", "GAM only high", "MaxEnt only high", "Both high")
)

rank_diff_r <- grid_template
values(rank_diff_r)[comparison_df$cell_id] <- comparison_df$rank_diff

writeRaster(agreement_r,   file.path(out_dir, "gam_maxent_agreement_top25.tif"), overwrite = TRUE)
writeRaster(gam_rank_r,    file.path(out_dir, "gam_rank_relative.tif"), overwrite = TRUE)
writeRaster(maxent_rank_r, file.path(out_dir, "maxent_rank_relative.tif"), overwrite = TRUE)
writeRaster(rank_diff_r,   file.path(out_dir, "gam_maxent_rank_diff.tif"), overwrite = TRUE)
message("✓ Derived rasters saved.")

# -------------------------------------------------------------------------
# STEP 11: Global correlation and multi-threshold concordance metrics
# -------------------------------------------------------------------------
message("--> Step 11: Calculating global statistics...")

cohen_kappa_2x2 <- function(gam_high, maxent_high) {
  tab <- table(gam_high, maxent_high)
  if (!all(dim(tab) == c(2, 2))) {
    full <- matrix(0, 2, 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
    full[rownames(tab), colnames(tab)] <- tab
    tab <- full
  }
  n  <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) / n * colSums(tab) / n)
  if ((1 - pe) == 0) return(NA_real_)
  (po - pe) / (1 - pe)
}

threshold_concordance <- function(df, pct) {
  g_thr <- quantile(df$gam, pct, na.rm = TRUE)
  m_thr <- quantile(df$maxent, pct, na.rm = TRUE)
  
  g_h <- df$gam > g_thr
  m_h <- df$maxent > m_thr
  
  both_h  <- g_h & m_h
  inter_n <- sum(both_h)
  union_n <- sum(g_h | m_h)
  
  jaccard  <- if (union_n > 0) inter_n / union_n else NA_real_
  sorensen <- if ((sum(g_h) + sum(m_h)) > 0) 2 * inter_n / (sum(g_h) + sum(m_h)) else NA_real_
  kappa    <- cohen_kappa_2x2(g_h, m_h)
  
  data.frame(
    threshold       = sprintf("Top %d%%", round((1 - pct) * 100)),
    n_pixels        = nrow(df),
    jaccard         = round(jaccard, 3),
    sorensen        = round(sorensen, 3),
    kappa           = round(kappa, 3),
    pct_both_high   = round(100 * mean(both_h), 1),
    pct_gam_only    = round(100 * mean(g_h & !m_h), 1),
    pct_maxent_only = round(100 * mean(!g_h & m_h), 1)
  )
}

spearman_res <- cor.test(comparison_df$gam, comparison_df$maxent, method = "spearman")
pearson_res  <- cor.test(comparison_df$gam, comparison_df$maxent, method = "pearson")

message(sprintf(
  "Spatial correlation (n = %s): Spearman rho = %.3f (p = %.4g); Pearson r = %.3f (p = %.4g)",
  format(nrow(comparison_df), big.mark = ","),
  unname(spearman_res$estimate), spearman_res$p.value,
  unname(pearson_res$estimate),  pearson_res$p.value
))

concordance_table <- bind_rows(
  threshold_concordance(comparison_df, 0.90),
  threshold_concordance(comparison_df, 0.75),
  threshold_concordance(comparison_df, 0.50)
)

print(concordance_table)
write.csv(
  concordance_table,
  file.path(out_dir, "gam_maxent_concordance_table.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# STEP 12: Publication-ready plotting theme
# -------------------------------------------------------------------------
theme_thesis_map <- function(base_family = "sans", base_size = 10) {
  theme_minimal(base_family = base_family, base_size = base_size) %+replace%
    theme(
      plot.title        = element_text(size = rel(1.10), face = "bold", color = "grey15", hjust = 0, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = rel(0.88), color = "grey35", hjust = 0, margin = margin(b = 5)),
      axis.title        = element_text(size = rel(0.82), face = "bold", color = "grey25"),
      axis.text         = element_text(size = rel(0.72), color = "grey35"),
      panel.grid        = element_line(color = "grey92", linewidth = 0.2),
      legend.title      = element_text(size = rel(0.78), face = "bold", color = "grey20"),
      legend.text       = element_text(size = rel(0.70), color = "grey30"),
      legend.position   = "right",
      legend.key.height = unit(0.42, "cm"),
      legend.key.width  = unit(0.32, "cm"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin       = margin(t = 5, r = 5, b = 5, l = 5)
    )
}

# -------------------------------------------------------------------------
# STEP 13: Build 4-panel figure with colour-blind friendly palettes
# -------------------------------------------------------------------------
message("--> Step 13: Rendering 4-panel figure...")

agreement_cols <- c(
  "Both low"        = "#D9D9D9",
  "GAM only high"   = "#0072B2",
  "MaxEnt only high"= "#E69F00",
  "Both high"       = "#009E73"
)

j_val <- concordance_table$jaccard[concordance_table$threshold == "Top 25%"]

p_gam <- ggplot() +
  geom_spatraster(data = gam_rank_r) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    name = "Relative\nsuitability",
    na.value = "transparent"
  ) +
  labs(
    title = "A. GAM relative suitability",
    subtitle = "Pakistan-calibrated"
  ) +
  guides(fill = guide_colorbar(barheight = unit(2.6, "cm"), barwidth = unit(0.38, "cm"))) +
  theme_thesis_map()

p_maxent <- ggplot() +
  geom_spatraster(data = maxent_rank_r) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    name = "Relative\nsuitability",
    na.value = "transparent"
  ) +
  labs(
    title = "B. MaxEnt relative suitability",
    subtitle = "India-calibrated"
  ) +
  guides(fill = guide_colorbar(barheight = unit(2.6, "cm"), barwidth = unit(0.38, "cm"))) +
  theme_thesis_map()

p_agreement <- ggplot() +
  geom_spatraster(data = agreement_factor) +
  scale_fill_manual(
    values = agreement_cols,
    name = "Concordance\nclass",
    na.translate = FALSE
  ) +
  labs(
    title = "C. Categorical spatial concordance",
    subtitle = sprintf("Top 25%% (J = %.3f)", j_val)
  ) +
  guides(fill = guide_legend(keyheight = unit(0.38, "cm"), keywidth = unit(0.38, "cm"))) +
  theme_thesis_map()

p_rank_diff <- ggplot() +
  geom_spatraster(data = rank_diff_r) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-1, 1),
    name = "GAM rank −\nMaxEnt rank",
    na.value = "transparent"
  ) +
  labs(
    title = "D. Relative rank difference",
    subtitle = "GAM percentile − MaxEnt percentile"
  ) +
  guides(fill = guide_colorbar(barheight = unit(2.6, "cm"), barwidth = unit(0.38, "cm"))) +
  theme_thesis_map()

final_figure <- (p_gam | p_maxent) / (p_agreement | p_rank_diff) +
  plot_annotation(
    title = "Figure 6. GAM–MaxEnt spatial concordance across Pakistan",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", color = "grey15")
    )
  )

ggsave(
  file.path(out_dir, "fig6_pakistan_concordance.png"),
  final_figure,
  width = 11,
  height = 9.5,
  dpi = 300,
  bg = "white"
)

message("✓ Final figure saved: ", file.path(out_dir, "gam_maxent_concordance_4panel.png"))
message("✓ Concordance table saved: ", file.path(out_dir, "gam_maxent_concordance_table.csv"))
message("✓ Improved workflow complete.")