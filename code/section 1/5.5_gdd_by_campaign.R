# ==============================================================================
# FIGURE 3: GDD RESPONSE CURVES + HABITAT RUG + CAMPAIGN COLOUR STRIPS
# ==============================================================================

library(ggplot2)
library(dplyr)
library(lubridate)
library(gratia)
library(sf)
library(scales)

# ------------------------------------------------------------------------------
# 0. Load
# ------------------------------------------------------------------------------
out_dir       <- "~/Documents/DISSERT/EDITS/section 1/results"
survey_data   <- read.csv(
  "~/Documents/DISSERT/EDITS/section 1/survey_data.csv",
  stringsAsFactors = FALSE
)
gam_gdd_k_200 <- readRDS(
  "~/Documents/DISSERT/EDITS/section 1/gam_rds/gam_gdd_k_200.rds"
)

# ------------------------------------------------------------------------------
# 1. Parse & prepare
# ------------------------------------------------------------------------------
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

# Projected coords
if (!all(c("x_km", "y_km") %in% colnames(survey_data))) {
  pts_proj         <- st_as_sf(survey_data,
                               coords = c("longitude", "latitude"),
                               crs    = 4326) %>%
    st_transform(32643)
  coords_m         <- st_coordinates(pts_proj)
  survey_data$x_km <- coords_m[, 1] / 1000
  survey_data$y_km <- coords_m[, 2] / 1000
}

# ------------------------------------------------------------------------------
# 2. Habitat strip labels with n
# ------------------------------------------------------------------------------
hab_order <- c(
  "Agricultural_disturbed",
  "Agricultural_managed",
  "Semi_natural_open",
  "Urban_disturbed",
  "Wet_habitat",
  "Woody"
)

hab_counts <- survey_data %>%
  filter(!is.na(hab_group)) %>%
  count(hab_group) %>%
  mutate(
    hab_group   = factor(hab_group, levels = hab_order),
    strip_label = paste0(gsub("_", " ", hab_group),
                         "\n(n = ", comma(n), ")")
  ) %>%
  arrange(hab_group)

# Ordered strip labels
strip_levels <- hab_counts$strip_label

survey_data <- survey_data %>%
  left_join(hab_counts %>% select(hab_group, strip_label), by = "hab_group") %>%
  mutate(strip_label = factor(strip_label, levels = strip_levels))

# ------------------------------------------------------------------------------
# 3. GAM prediction grid
# ------------------------------------------------------------------------------
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
) %>%
  left_join(hab_counts %>% select(hab_group, strip_label), by = "hab_group") %>%
  mutate(strip_label = factor(strip_label, levels = strip_levels))

# ------------------------------------------------------------------------------
# 4. Campaign colours
# ------------------------------------------------------------------------------
campaign_colors <- c(
  "C1: Dec–Jan" = "#E69F00",
  "C2: Feb–Mar" = "#56B4E9",
  "C3: Jun–Aug" = "#009E73",
  "C4: Oct–Dec" = "#CC79A7"
)

# ------------------------------------------------------------------------------
# 5. Rug data subsets
# ------------------------------------------------------------------------------
x_limits <- range(survey_data$gdd_30d, na.rm = TRUE)
x_expand <- expansion(mult = 0.01)

rug_all <- survey_data %>%
  filter(!is.na(gdd_30d), !is.na(strip_label))

rug_c1 <- rug_all %>% filter(campaign == "C1: Dec–Jan")
rug_c2 <- rug_all %>% filter(campaign == "C2: Feb–Mar")
rug_c3 <- rug_all %>% filter(campaign == "C3: Jun–Aug")
rug_c4 <- rug_all %>% filter(campaign == "C4: Oct–Dec")

# ------------------------------------------------------------------------------
# 6. Dummy data for legend — proper x and y supplied
# ------------------------------------------------------------------------------
dummy_legend <- data.frame(
  gdd_30d    = rep(mean(x_limits), 4),   # valid x value
  .fitted    = rep(0, 4),                # valid y value
  strip_label = factor(strip_levels[1], levels = strip_levels),
  campaign    = factor(names(campaign_colors),
                       levels = names(campaign_colors))
)

# ------------------------------------------------------------------------------
# 7. Build figure
# ------------------------------------------------------------------------------
fig_final <- ggplot(
  gdd_hab_fv,
  aes(x = gdd_30d, y = .fitted)
) +
  
  # CI ribbon
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha  = 0.15,
    fill   = "#2166AC",
    colour = NA
  ) +
  
  # Fitted curve
  geom_line(
    colour    = "#2166AC",
    linewidth = 1.0
  ) +
  
  # Black rug — all observations
  geom_rug(
    data        = rug_all,
    aes(x       = gdd_30d),
    sides       = "b",
    alpha       = 0.3,
    linewidth   = 0.3,
    colour      = "black",
    length      = unit(0.04, "npc"),
    inherit.aes = FALSE
  ) +
  
  # Campaign 1 rug
  geom_rug(
    data        = rug_c1,
    aes(x       = gdd_30d),
    sides       = "b",
    colour      = "#E69F00",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    inherit.aes = FALSE
  ) +
  
  # Campaign 2 rug
  geom_rug(
    data        = rug_c2,
    aes(x       = gdd_30d),
    sides       = "b",
    colour      = "#56B4E9",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    inherit.aes = FALSE
  ) +
  
  # Campaign 3 rug
  geom_rug(
    data        = rug_c3,
    aes(x       = gdd_30d),
    sides       = "b",
    colour      = "#009E73",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    inherit.aes = FALSE
  ) +
  
  # Campaign 4 rug
  geom_rug(
    data        = rug_c4,
    aes(x       = gdd_30d),
    sides       = "b",
    colour      = "#CC79A7",
    alpha       = 0.7,
    linewidth   = 0.4,
    length      = unit(0.04, "npc"),
    inherit.aes = FALSE
  ) +
  
  # Dummy layer for legend only — invisible points
  geom_point(
    data        = dummy_legend,
    aes(x       = gdd_30d,
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
      override.aes = list(
        size      = 4,
        alpha     = 1,
        shape     = 15    # filled square — clear colour swatch
      )
    )
  ) +
  
  # Facet
  facet_wrap(~ strip_label, ncol = 3) +
  
  # Axes
  scale_x_continuous(limits = x_limits, expand = x_expand) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = percent_format(accuracy = 1)
  ) +
  
  labs(
    title = "Figure 3: GDD response by habitat",
    x     = "30-day growing degree days (GDD)",
    y     = "Predicted probability"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    axis.title      = element_text(face = "bold", size = 10),
    axis.text       = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    strip.text      = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 9),
    plot.margin     = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------------------------
ggsave(
  file.path(out_dir, "fig3_hab_campaign_rug.png"),
  fig_final,
  width  = 11,
  height = 8,
  dpi    = 300
)

cat("✅ Saved: fig3_hab_campaign_rug.png\n")