# ==============================================================================
# SCRIPT 05: FINAL PUBLICATION FIGURES & RESIDUAL MAPS (REVIEWER-ALIGNED)
# ==============================================================================

library(tidyverse)
library(mgcv)
library(gratia)
library(marginaleffects)
library(patchwork)
library(gstat)
library(sf)
library(grid)

out_dir <- "/home/anaga-ambady/Documents/DISSERT/EDITS/section 1/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 0. LOAD CLEAN DATA + MODEL
# ------------------------------------------------------------------------------

survey_data <- read.csv(
  "/home/anaga-ambady/Documents/DISSERT/DATA/MODEL/results/survey_data.csv"
) %>%
  dplyr::mutate(
    hab_group = factor(
      hab_group,
      levels = c(
        "Agricultural_disturbed",
        "Agricultural_managed",
        "Semi_natural_open",
        "Urban_disturbed",
        "Wet_habitat",
        "Woody"
      )
    ),
    road_type = factor(road_type),
    plotID    = factor(plotID)
  )

gam_gdd_k_200 <- readRDS(
  "/home/anaga-ambady/Documents/DISSERT/EDITS/section 1/gam_rds/gam_gdd_k_200.rds"
)

# --- Ensure projected coordinates exist and match model CRS (EPSG:32642) ---
if (!all(c("x_km", "y_km") %in% colnames(survey_data))) {
  pts_proj <- sf::st_as_sf(survey_data, coords = c("longitude", "latitude"), crs = 4326) %>%
    sf::st_transform(32642)
  coords_m <- sf::st_coordinates(pts_proj)
  survey_data$x_km <- coords_m[, 1] / 1000
  survey_data$y_km <- coords_m[, 2] / 1000
}

# Drop rows with missing model variables to keep predictions and residuals aligned
model_vars <- c(
  "presence", "hab_group", "road_type", "precip_sum_30d",
  "gdd_30d", "altitude", "hab_perc", "x_km", "y_km"
)

survey_data <- survey_data %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(model_vars), ~ !is.na(.)))

# ------------------------------------------------------------------------------
# Helpers: build complete prediction data for this GAM
# ------------------------------------------------------------------------------

get_ref_road <- function(survey_df) {
  if ("District Road" %in% levels(survey_df$road_type)) {
    "District Road"
  } else {
    levels(survey_df$road_type)[1]
  }
}

build_pred_grid <- function(gdd_values, hab_levels, survey_df) {
  ref_road <- get_ref_road(survey_df)
  
  tidyr::expand_grid(
    gdd_30d   = gdd_values,
    hab_group = hab_levels
  ) %>%
    dplyr::mutate(
      precip_sum_30d = median(survey_df$precip_sum_30d, na.rm = TRUE),
      altitude       = median(survey_df$altitude, na.rm = TRUE),
      hab_perc       = median(survey_df$hab_perc, na.rm = TRUE),
      road_type      = factor(ref_road, levels = levels(survey_df$road_type)),
      x_km           = median(survey_df$x_km, na.rm = TRUE),
      y_km           = median(survey_df$y_km, na.rm = TRUE)
    )
}

build_single_predict_df <- function(x_seq, focal_var, survey_df) {
  ref_road <- get_ref_road(survey_df)
  ref_hab  <- levels(survey_df$hab_group)[1]
  
  out <- tibble::tibble(x_tmp = x_seq)
  names(out)[1] <- focal_var
  
  out %>%
    dplyr::mutate(
      hab_group      = factor(ref_hab, levels = levels(survey_df$hab_group)),
      road_type      = factor(ref_road, levels = levels(survey_df$road_type)),
      precip_sum_30d = median(survey_df$precip_sum_30d, na.rm = TRUE),
      gdd_30d        = median(survey_df$gdd_30d, na.rm = TRUE),
      altitude       = median(survey_df$altitude, na.rm = TRUE),
      hab_perc       = median(survey_df$hab_perc, na.rm = TRUE),
      x_km           = median(survey_df$x_km, na.rm = TRUE),
      y_km           = median(survey_df$y_km, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      precip_sum_30d = ifelse(focal_var == "precip_sum_30d", .data[[focal_var]], precip_sum_30d),
      gdd_30d        = ifelse(focal_var == "gdd_30d", .data[[focal_var]], gdd_30d)
    )
}

# ==============================================================================
# 1. PARAMETRIC FOREST PLOTS & TABULAR SUMMARY
# ==============================================================================

plot_df <- data.frame(
  Term = c(
    "ROAD TYPE",
    "Paths / Trails", 
    "Not Paved", 
    "Local Roads", 
    "National Highways", 
    "HABITAT GROUP",
    "Agricultural Managed", 
    "Woody", 
    "Urban Disturbed", 
    "Wet Habitat", 
    "Semi-natural Open"
  ),
  Estimate = c(NA, -1.0766945, -0.6124466, -0.3395976, 0.0394727, NA, -1.5686065, -0.7720707, -0.5991479, -0.5796283, 0.2374386),
  SE       = c(NA, 0.2554015, 0.1675768, 0.1325969, 0.1673686, NA, 0.1610074, 0.1742162, 0.2953930, 0.2335589, 0.1375685),
  p_val    = c(NA, 2.49e-05, 2.57e-04, 0.0104, 0.8136, NA, 1.98e-22, 9.35e-06, 0.0425, 0.0131, 0.0844),
  is_header = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
) %>%
  dplyr::mutate(
    Term         = factor(Term, levels = Term),
    low          = Estimate - 1.96 * SE,
    high         = Estimate + 1.96 * SE,
    OR           = exp(Estimate),
    OR_low       = exp(low),
    OR_high      = exp(high),
    OR_CI_string = dplyr::if_else(is_header, "", sprintf("%.2f [%.2f, %.2f]", OR, OR_low, OR_high)),
    p_string     = dplyr::if_else(is_header, "", dplyr::if_else(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val))),
    Significance = dplyr::case_when(
      is_header ~ NA_character_,
      p_val < 0.05 ~ "Significant (p < 0.05)",
      TRUE ~ "Non-Significant (p ≥ 0.05)"
    )
  )

y_fontfaces <- dplyr::if_else(levels(plot_df$Term) %in% c("HABITAT GROUP", "ROAD TYPE"), "bold", "plain")

# Forest plot panel with expanded y-axis padding
p_forest <- ggplot(plot_df, aes(x = Estimate, y = Term, xmin = low, xmax = high, color = Significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.8, na.rm = TRUE) +
  scale_color_manual(
    values = c("Significant (p < 0.05)" = "#0072B2", "Non-Significant (p ≥ 0.05)" = "#999999"),
    na.value = "transparent"
  ) +
  scale_x_continuous(limits = c(-2.2, 0.8), breaks = seq(-2, 0.5, by = 0.5)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) + # Adds proportional padding above/below rows
  labs(
    x = "Log-Odds Coefficient (95% CI)",
    y = NULL,
    title = "Figure 2: Parametric Predictor Effects on Parthenium Presence",
    subtitle = "Baselines: Agricultural Disturbed (Habitat) | District Road (Road Type)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey30", margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(face = y_fontfaces, size = 10.5, color = "black"),
    legend.position  = "bottom",
    legend.title     = element_blank()
  )

# Tabular summary panel aligned with forest y-axis
p_table <- ggplot(plot_df, aes(y = Term)) +
  geom_text(aes(x = 0, label = OR_CI_string), hjust = 0.5, size = 3.6) +
  geom_text(aes(x = 1, label = p_string), hjust = 0.5, size = 3.6) +
  annotate("text", x = 0, y = 0.1, label = "Odds Ratio [95% CI]", fontface = "bold", size = 3.8, hjust = 0.5) +
  annotate("text", x = 1, y = 0.1, label = "p-value", fontface = "bold", size = 3.8, hjust = 0.5) +
  scale_x_continuous(limits = c(-0.5, 1.5)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) + # Matches forest panel padding
  coord_cartesian(ylim = c(1, 11), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(t = 42, r = 10, b = 38, l = 10)
  )

combined_figure <- p_forest + p_table + plot_layout(widths = c(2.3, 1))

# Save with height = 7.5 to increase row spacing across all 11 rows
ggsave(file.path(out_dir, "fig2_parametric_forest_table.png"), combined_figure, width = 11, height = 7.5, dpi = 300)

# ==============================================================================
# 2. PRECIPITATION RESPONSE CURVE
# ============================================================================

precip_seq <- seq(
  min(survey_data$precip_sum_30d, na.rm = TRUE),
  max(survey_data$precip_sum_30d, na.rm = TRUE),
  length.out = 200
)

precip_slice <- build_single_predict_df(
  x_seq = precip_seq,
  focal_var = "precip_sum_30d",
  survey_df = survey_data
)

precip_fv <- fitted_values(
  gam_gdd_k_200,
  data = precip_slice,
  scale = "response",
  exclude = "s(x_km,y_km)"
)

fig3_precip <- ggplot(precip_fv, aes(x = precip_sum_30d, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#009E73", alpha = 0.16) +
  geom_line(colour = "#00664A", linewidth = 1.2) +
  geom_rug(
    data = survey_data,
    aes(x = precip_sum_30d),
    inherit.aes = FALSE,
    colour = "#00664A",
    alpha = 0.18,
    sides = "b"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x = "30-day cumulative precipitation (mm)",
    y = "Predicted probability",
    title = "Figure 4: Response of Parthenium to 30-day precipitation",
    subtitle = "Global GAM curve with 95% intervals; other predictors held at typical values"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "fig3_precipitation_smooth.png"), fig3_precip, width = 9, height = 6, dpi = 300)

# ==============================================================================
# 3. POPULATION-AVERAGED GDD RESPONSE CURVE
# ============================================================================

gdd_seq_global <- seq(
  min(survey_data$gdd_30d, na.rm = TRUE),
  max(survey_data$gdd_30d, na.rm = TRUE),
  length.out = 200
)

gdd_slice <- build_single_predict_df(
  x_seq = gdd_seq_global,
  focal_var = "gdd_30d",
  survey_df = survey_data
)

gdd_fv <- fitted_values(
  gam_gdd_k_200,
  data = gdd_slice,
  scale = "response",
  exclude = "s(x_km,y_km)"
)

fig4_gdd_pop <- ggplot(gdd_fv, aes(x = gdd_30d, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#0072B2", alpha = 0.16) +
  geom_line(colour = "#004B75", linewidth = 1.2) +
  geom_rug(
    data = survey_data,
    aes(x = gdd_30d),
    inherit.aes = FALSE,
    colour = "#004B75",
    alpha = 0.18,
    sides = "b"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x = "30-day growing degree days (GDD)",
    y = "Predicted probability",
    title = "Figure 4: Global response of Parthenium to growing degree days",
    subtitle = "Global GAM curve with 95% intervals; other predictors held at typical values"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "fig4_gdd_population_smooth.png"), fig4_gdd_pop, width = 9, height = 6, dpi = 300)

# ==============================================================================
# 4. BY-HABITAT GDD RESPONSE CURVES (REVIEWER-ALIGNED)
# ============================================================================

hab_labels <- c(
  "Agricultural_disturbed" = "Agricultural disturbed",
  "Agricultural_managed"   = "Agricultural managed",
  "Semi_natural_open"      = "Semi-natural open",
  "Urban_disturbed"        = "Urban disturbed",
  "Wet_habitat"            = "Wet habitat",
  "Woody"                  = "Woody"
)

# Summarise GDD support by habitat
hab_support <- survey_data %>%
  dplyr::filter(!is.na(gdd_30d), !is.na(hab_group)) %>%
  dplyr::group_by(hab_group) %>%
  dplyr::summarise(
    gdd_min = min(gdd_30d, na.rm = TRUE),
    gdd_max = max(gdd_30d, na.rm = TRUE),
    gdd_q05 = quantile(gdd_30d, 0.05, na.rm = TRUE),
    gdd_q95 = quantile(gdd_30d, 0.95, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

# Add sample sizes to facet labels
hab_labels_n <- setNames(
  paste0(hab_labels[as.character(hab_support$hab_group)], "\n(n = ", hab_support$n, ")"),
  as.character(hab_support$hab_group)
)

# Build habitat-specific prediction grid ONLY for the central supported range
ref_road <- get_ref_road(survey_data)
ref_precip <- median(survey_data$precip_sum_30d, na.rm = TRUE)
ref_alt    <- median(survey_data$altitude, na.rm = TRUE)
ref_habpct <- median(survey_data$hab_perc, na.rm = TRUE)
ref_x      <- median(survey_data$x_km, na.rm = TRUE)
ref_y      <- median(survey_data$y_km, na.rm = TRUE)

gdd_hab_slice <- purrr::map_dfr(seq_len(nrow(hab_support)), function(i) {
  this_hab <- hab_support$hab_group[i]
  this_q05 <- hab_support$gdd_q05[i]
  this_q95 <- hab_support$gdd_q95[i]
  
  tibble::tibble(
    gdd_30d = seq(this_q05, this_q95, length.out = 200),
    hab_group = factor(this_hab, levels = levels(survey_data$hab_group)),
    precip_sum_30d = ref_precip,
    altitude = ref_alt,
    hab_perc = ref_habpct,
    road_type = factor(ref_road, levels = levels(survey_data$road_type)),
    x_km = ref_x,
    y_km = ref_y
  )
})

gdd_hab_fv <- fitted_values(
  gam_gdd_k_200,
  data = gdd_hab_slice,
  scale = "response",
  exclude = "s(x_km,y_km)"
) %>%
  dplyr::left_join(hab_support, by = "hab_group")

fig2_gdd_by_habitat <- ggplot() +
  geom_rect(
    data = hab_support,
    aes(xmin = gdd_q05, xmax = gdd_q95, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey96",
    alpha = 1
  ) +
  geom_line(
    data = gdd_hab_fv,
    aes(x = gdd_30d, y = .fitted, group = hab_group),
    colour = "#0072B2",
    linewidth = 0.9
  ) +
  geom_rug(
    data = survey_data %>% dplyr::filter(!is.na(gdd_30d), !is.na(hab_group)),
    aes(x = gdd_30d),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.10,
    colour = "black"
  ) +
  facet_wrap(~ hab_group, ncol = 3, labeller = as_labeller(hab_labels_n)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 0.40)) +
  labs(
    x = "30-day growing degree days (GDD)",
    y = "Predicted probability",
    title = "Figure 3: GDD response by habitat",
    subtitle = "Blue lines show habitat-specific fitted responses only within the central observed GDD range; rugs show sample support"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 10.5, color = "grey30", margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )

print(fig2_gdd_by_habitat)

ggsave(
  file.path(out_dir, "fig2_gdd_by_habitat_reviewer_aligned.png"),
  fig2_gdd_by_habitat,
  width = 11,
  height = 7.5,
  dpi = 300
)

# ==============================================================================
# 5. SPATIAL RESIDUAL MAP & VARIOGRAM
# ============================================================================

survey_data$resid_dev <- residuals(gam_gdd_k_200, type = "deviance")

fig_map <- ggplot(survey_data, aes(x = x_km, y = y_km, colour = resid_dev)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_viridis_c(option = "plasma", name = "Deviance\nresidual") +
  theme_minimal(base_size = 12) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Spatial map of model deviance residuals",
    x = "Easting (km)",
    y = "Northing (km)"
  )

ggsave(file.path(out_dir, "spatial_residual_map.png"), fig_map, width = 8, height = 6, dpi = 300)

pts_sp <- survey_data %>%
  dplyr::mutate(
    pearson_resid = as.numeric(residuals(gam_gdd_k_200, type = "pearson"))
  ) %>%
  dplyr::group_by(x_km, y_km) %>%
  dplyr::summarise(
    resid = mean(pearson_resid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  sf::st_as_sf(coords = c("x_km", "y_km"), crs = 32642) %>%
  as("Spatial")

vg <- variogram(resid ~ 1, pts_sp)
png(file.path(out_dir, "variogram_residuals.png"), width = 800, height = 600)
plot(vg, main = "Variogram of aggregated Pearson residuals", xlab = "Distance (km)", ylab = "Semivariance")
dev.off()
