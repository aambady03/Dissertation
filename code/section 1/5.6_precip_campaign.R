# ==============================================================================
# SCRIPT 05: FINAL PUBLICATION FIGURES & RESIDUAL MAPS
# ==============================================================================

library(tidyverse)
library(mgcv)
library(gratia)
library(marginaleffects)
library(patchwork)
library(gstat)
library(sf)

out_dir       <- "~/Documents/DISSERT/EDITS/section 1/results"
survey_data   <- read.csv("~/Documents/DISSERT/EDITS/section 1/survey_data.csv")
gam_gdd_k_200 <- readRDS("~/Documents/DISSERT/EDITS/section 1/gam_rds/gam_gdd_k_200.rds")

# --- 0. Ensure Projected Coordinates ---
if (!all(c("x_km", "y_km") %in% colnames(survey_data))) {
  pts_proj <- st_as_sf(survey_data, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(32643)
  coords_m <- st_coordinates(pts_proj)
  survey_data$x_km <- coords_m[, 1] / 1000
  survey_data$y_km <- coords_m[, 2] / 1000
}

# --- 0b. Parse dates & campaign labels (used in Section 2) ---
survey_data <- survey_data %>%
  mutate(
    date_parsed = as.Date(today, format = "%m/%d/%Y %H:%M"),
    campaign    = factor(
      campaign,
      levels = sort(unique(campaign)),
      labels = c(
        "C1: Dec–Jan",
        "C2: Feb–Mar",
        "C3: Jun–Aug",
        "C4: Oct–Dec"
      )
    )
  )

campaign_colors <- c(
  "C1: Dec–Jan" = "#E69F00",
  "C2: Feb–Mar" = "#56B4E9",
  "C3: Jun–Aug" = "#009E73",
  "C4: Oct–Dec" = "#CC79A7"
)

# ==============================================================================
# 1. PARAMETRIC FOREST PLOTS & TABULAR SUMMARY
# ==============================================================================

plot_df <- data.frame(
  Term = c(
    "Agricultural Managed", "Semi-natural Open", "Urban Disturbed", "Wet Habitat", "Woody",
    "Local Roads", "National Highways", "Not Paved", "Paths / Trails"
  ),
  Estimate = c(-1.5686065, 0.2374386, -0.5991479, -0.5796283, -0.7720707,
               -0.3395976,  0.0394727, -0.6124466, -1.0766945),
  SE       = c( 0.1610074, 0.1375685,  0.2953930,  0.2335589,  0.1742162,
                0.1325969,  0.1673686,  0.1675768,  0.2554015),
  p_val    = c(1.98e-22, 0.0844, 0.0425, 0.0131, 9.35e-06,
               0.0104,   0.8136, 2.57e-04, 2.49e-05)
) %>%
  mutate(
    low          = Estimate - 1.96 * SE,
    high         = Estimate + 1.96 * SE,
    OR           = exp(Estimate),
    OR_low       = exp(low),
    OR_high      = exp(high),
    OR_CI_string = sprintf("%.2f [%.2f, %.2f]", OR, OR_low, OR_high),
    p_string     = ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    Significance = factor(if_else(p_val < 0.05,
                                  "Significant (p < 0.05)",
                                  "Non-Significant (p ≥ 0.05)")),
    Term         = reorder(Term, Estimate)
  )

p_forest <- ggplot(
  plot_df,
  aes(x = Estimate, y = Term, xmin = low, xmax = high, color = Significance)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(
    "Significant (p < 0.05)"     = "#0072B2",
    "Non-Significant (p ≥ 0.05)" = "#999999"
  )) +
  scale_x_continuous(limits = c(-2.0, 0.8), breaks = seq(-2, 0.5, by = 0.5)) +
  labs(
    x        = "Log-Odds Coefficient (95% CI)",
    y        = NULL,
    title    = "Figure 1: Parametric Predictor Effects on Parthenium Presence",
    subtitle = "Baselines: Agricultural Disturbed (Habitat) | District Road (Road Type)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey30", margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 10.5, color = "black"),
    legend.position  = "bottom",
    legend.title     = element_blank()
  )

p_table <- ggplot(plot_df, aes(y = Term)) +
  geom_text(aes(x = 0, label = OR_CI_string), hjust = 0.5, size = 3.6) +
  geom_text(aes(x = 1, label = p_string),     hjust = 0.5, size = 3.6) +
  scale_x_continuous(
    limits = c(-0.5, 1.5),
    breaks = c(0, 1),
    labels = c("Odds Ratio [95% CI]", "p-value")
  ) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10.5, color = "black"),
    plot.margin = margin(t = 38, r = 10, b = 40, l = 10)
  )

combined_figure <- p_forest + p_table + plot_layout(widths = c(2.3, 1))
ggsave(
  file.path(out_dir, "fig1_parametric_forest_table.png"),
  combined_figure, width = 11, height = 6, dpi = 300
)

# ==============================================================================
# 2. PRECIPITATION RESPONSE CURVE + CAMPAIGN RUGS
# ==============================================================================

# --- Prediction grid via data_slice (same as your original) ---
precip_slice <- data_slice(
  gam_gdd_k_200,
  precip_sum_30d = evenly(survey_data$precip_sum_30d, n = 200)
)

precip_fv <- fitted_values(
  gam_gdd_k_200,
  data    = precip_slice,
  scale   = "response",
  exclude = "s(x_km,y_km)"
)

# --- Rug subsets by campaign ---
rug_base <- survey_data %>% filter(!is.na(precip_sum_30d))
rug_c1   <- rug_base %>% filter(campaign == "C1: Dec–Jan")
rug_c2   <- rug_base %>% filter(campaign == "C2: Feb–Mar")
rug_c3   <- rug_base %>% filter(campaign == "C3: Jun–Aug")
rug_c4   <- rug_base %>% filter(campaign == "C4: Oct–Dec")

# --- Dummy data for legend ---
dummy_precip <- data.frame(
  precip_sum_30d = rep(median(survey_data$precip_sum_30d, na.rm = TRUE), 4),
  .fitted        = rep(0, 4),
  campaign       = factor(names(campaign_colors), levels = names(campaign_colors))
)

# --- Plot ---
fig2_precip <- ggplot(
  precip_fv,
  aes(x = precip_sum_30d, y = .fitted)
) +
  
  # CI ribbon
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),
    fill   = "#009E73",
    alpha  = 0.2,
    colour = NA
  ) +
  
  # Fitted curve
  geom_line(colour = "#00664A", linewidth = 1.3) +
  
  # Black rug — all observations
  geom_rug(
    data        = rug_base,
    aes(x       = precip_sum_30d),
    inherit.aes = FALSE,
    colour      = "black",
    alpha       = 0.18,
    linewidth   = 0.3,
    length      = unit(0.04, "npc"),
    sides       = "b"
  ) +
  
  # Campaign 1
  geom_rug(
    data        = rug_c1,
    aes(x       = precip_sum_30d),
    inherit.aes = FALSE,
    colour      = "#E69F00",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    sides       = "b"
  ) +
  
  # Campaign 2
  geom_rug(
    data        = rug_c2,
    aes(x       = precip_sum_30d),
    inherit.aes = FALSE,
    colour      = "#56B4E9",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    sides       = "b"
  ) +
  
  # Campaign 3
  geom_rug(
    data        = rug_c3,
    aes(x       = precip_sum_30d),
    inherit.aes = FALSE,
    colour      = "#009E73",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    sides       = "b"
  ) +
  
  # Campaign 4
  geom_rug(
    data        = rug_c4,
    aes(x       = precip_sum_30d),
    inherit.aes = FALSE,
    colour      = "#CC79A7",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    sides       = "b"
  ) +
  
  # Dummy layer — legend only, invisible on plot
  geom_point(
    data        = dummy_precip,
    aes(x       = precip_sum_30d,
        y       = .fitted,
        colour  = campaign),
    size        = 0,
    alpha       = 0,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = campaign_colors,
    name   = "Campaign",
    guide  = guide_legend(
      override.aes = list(size = 4, alpha = 1, shape = 15)
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  
  labs(
    title = "Figure 4: Response of Parthenium to 30-day Precipitation",
    x     = "30-day Cumulative Precipitation (mm)",
    y     = "Predicted Probability"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = "black"),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9)
  )

ggsave(
  file.path(out_dir, "fig3_precipitation_smooth.png"),
  fig2_precip, width = 9, height = 6, dpi = 300
)

# ==============================================================================
# 3. POPULATION-AVERAGED GDD RESPONSE CURVE
# ==============================================================================

gdd_slice <- data_slice(gam_gdd_k_200, gdd_30d = evenly(survey_data$gdd_30d, n = 200))
gdd_fv    <- fitted_values(gam_gdd_k_200, data = gdd_slice, scale = "response", exclude = "s(x_km,y_km)")

fig3_gdd_pop <- ggplot(gdd_fv, aes(x = gdd_30d, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#0072B2", alpha = 0.2) +
  geom_line(colour = "#004B75", linewidth = 1.3) +
  geom_rug(
    data = survey_data, aes(x = gdd_30d),
    inherit.aes = FALSE, colour = "#004B75", alpha = 0.4, sides = "b"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "30-day Growing Degree Days (GDD)",
    y        = "Predicted Probability",
    title    = "Figure 3: Global Response of Parthenium to Growing Degree Days",
    subtitle = "Global GAM curve with 95% Bayesian credible intervals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "fig3_gdd_population_smooth.png"),
  fig3_gdd_pop, width = 9, height = 6, dpi = 300
)

# ==============================================================================
# 4. BY-HABITAT GDD RESPONSE CURVES
# ==============================================================================

okabe_ito_colors <- c(
  "Agricultural_disturbed" = "#E69F00",
  "Agricultural_managed"   = "#56B4E9",
  "Semi_natural_open"      = "#009E73",
  "Urban_disturbed"        = "#CC79A7",
  "Wet_habitat"            = "#D55E00",
  "Woody"                  = "#F0E442"
)

line_types <- c(
  "Agricultural_disturbed" = "solid",
  "Agricultural_managed"   = "solid",
  "Semi_natural_open"      = "solid",
  "Urban_disturbed"        = "solid",
  "Wet_habitat"            = "solid",
  "Woody"                  = "solid"
)

gdd_hab_slice <- data_slice(
  gam_gdd_k_200,
  gdd_30d   = evenly(survey_data$gdd_30d, n = 200),
  hab_group = unique(survey_data$hab_group)
)

gdd_hab_fv <- fitted_values(
  gam_gdd_k_200,
  data    = gdd_hab_slice,
  scale   = "response",
  exclude = "s(x_km,y_km)"
)

fig4_by_habitat <- ggplot(
  gdd_hab_fv,
  aes(x = gdd_30d, y = .fitted,
      colour = hab_group, fill = hab_group, linetype = hab_group)
) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = okabe_ito_colors, name = "Habitat Group") +
  scale_fill_manual(values = okabe_ito_colors,   name = "Habitat Group") +
  scale_linetype_manual(values = line_types,     name = "Habitat Group") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "30-day GDD",
    y        = "Predicted Probability",
    title    = "Figure 2: GDD Response Curves by Habitat Group",
    subtitle = "Fitted GAM curves evaluated across habitat interactions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave(
  file.path(out_dir, "fig2_gdd_by_habitat.png"),
  fig4_by_habitat, width = 9, height = 6, dpi = 300
)

# ==============================================================================
# 5. SPATIAL RESIDUAL MAP & VARIOGRAM
# ==============================================================================

survey_data$resid_dev <- residuals(gam_gdd_k_200, type = "deviance")

fig_map <- ggplot(survey_data, aes(x = x_km, y = y_km, colour = resid_dev)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_viridis_c(option = "plasma", name = "Deviance\nResidual") +
  theme_minimal(base_size = 12) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Spatial Map of Model Deviance Residuals",
    x     = "Easting (km)",
    y     = "Northing (km)"
  )

ggsave(
  file.path(out_dir, "spatial_residual_map.png"),
  fig_map, width = 8, height = 6, dpi = 300
)

pts_sp <- survey_data %>%
  dplyr::mutate(
    pearson_resid = as.numeric(residuals(gam_gdd_k_200, type = "pearson"))
  ) %>%
  dplyr::group_by(x_km, y_km) %>%
  dplyr::summarize(resid = mean(pearson_resid, na.rm = TRUE), .groups = "drop") %>%
  st_as_sf(coords = c("x_km", "y_km")) %>%
  as("Spatial")

vg <- variogram(resid ~ 1, pts_sp)
png(file.path(out_dir, "variogram_residuals.png"), width = 800, height = 600)
plot(vg,
     main  = "Variogram of Aggregated Pearson Residuals",
     xlab  = "Distance (km)",
     ylab  = "Semivariance")
dev.off()