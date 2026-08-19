library(ggplot2)
library(patchwork)
library(terra)
library(tidyterra)
library(grid)
library(cowplot)

# =========================================================================
# 1. Standardized Base Theme for 3x3 Spatial Grid (A4 Landscape)
# =========================================================================
theme_landscape_map <- function(base_size = 8.5) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title        = element_text(size = rel(1.05), face = "bold", color = "grey15", hjust = 0, margin = margin(b = 2)),
      plot.subtitle     = element_text(size = rel(0.85), color = "grey35", hjust = 0, margin = margin(b = 4)),
      axis.title        = element_blank(),
      axis.text         = element_text(size = rel(0.68), color = "grey40"),
      panel.grid        = element_line(color = "grey93", linewidth = 0.2),
      panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.35),
      legend.title      = element_text(size = rel(0.82), face = "bold", color = "grey20"),
      legend.text       = element_text(size = rel(0.72), color = "grey30"),
      legend.position   = "right",
      legend.key.height = unit(0.42, "cm"),
      legend.key.width  = unit(0.32, "cm"),
      legend.margin     = margin(l = 2, r = 2),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
    )
}

# =========================================================================
# 2. Plotting Functions with Uniform Legends
# =========================================================================
plot_suitability <- function(r, title_str, subtitle_str = NULL) {
  ggplot() +
    geom_spatraster(data = r) +
    scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      name = "Suitability\n(cloglog)",
      na.value = "transparent"
    ) +
    labs(title = title_str, subtitle = subtitle_str) +
    guides(fill = guide_colorbar(
      barheight = unit(2.8, "cm"),
      barwidth  = unit(0.32, "cm"),
      ticks.linewidth = 0.3,
      frame.colour = "grey50",
      frame.linewidth = 0.3
    )) +
    theme_landscape_map()
}

plot_overlap <- function(fact_r, title_str, subtitle_str = NULL) {
  overlap_cols <- c(
    "Neither suitable" = "#D9D9D9",
    "Parthenium only" = "#E69F00",
    "Beetle only"     = "#0072B2",
    "Both suitable"   = "#009E73"
  )
  
  ggplot() +
    geom_spatraster(data = fact_r) +
    scale_fill_manual(
      values = overlap_cols,
      name = "Categorical\nOverlap",
      na.translate = FALSE,
      drop = FALSE
    ) +
    labs(title = title_str, subtitle = subtitle_str) +
    guides(fill = guide_legend(
      keyheight = unit(0.36, "cm"),
      keywidth  = unit(0.36, "cm"),
      override.aes = list(color = "grey40", linewidth = 0.2)
    )) +
    theme_landscape_map()
}

# =========================================================================
# 3. Generate 3x3 Spatial Grid (Panels A–I)
# =========================================================================
pA <- plot_suitability(pair_now$parthenium, "A. Parthenium suitability", "Current")
pB <- plot_suitability(pair_30$parthenium,  "B. Parthenium suitability", "2030 (SSP2-4.5)")
pC <- plot_suitability(pair_50$parthenium,  "C. Parthenium suitability", "2050 (SSP2-4.5)")

pD <- plot_suitability(pair_now$beetle, "D. Beetle thermal proxy", "Current")
pE <- plot_suitability(pair_30$beetle,  "E. Beetle thermal proxy", "2030 (SSP2-4.5)")
pF <- plot_suitability(pair_50$beetle,  "F. Beetle thermal proxy", "2050 (SSP2-4.5)")

pG <- plot_overlap(pair_now$overlap_factor, "G. Niche overlap", "Current")
pH <- plot_overlap(pair_30$overlap_factor,  "H. Niche overlap", "2030 (SSP2-4.5)")
pI <- plot_overlap(pair_50$overlap_factor,  "I. Niche overlap", "2050 (SSP2-4.5)")

figure_3x3_clean <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)

# =========================================================================
# 4. Panel J: 2050 Summer Heat Stress with Fully Custom Legend
# =========================================================================
heat_cols <- c(
  "<= 38°C: within tolerance" = "#0072B2",
  "38-40°C: heat stress"      = "#E69F00",
  "> 40°C: likely unsuitable" = "#D55E00"
)

# Title block inside Panel J
pJ_title <- ggplot() +
  annotate(
    "text",
    x = 0, y = 1,
    label = "J. 2050 Summer Heat Stress",
    hjust = 0, vjust = 1,
    size = 4.0,
    fontface = "bold",
    color = "#8B0000"
  ) +
  annotate(
    "text",
    x = 0, y = 0.70,
    label = "Physiological upper limit filter for Z. bicolorata",
    hjust = 0, vjust = 1,
    size = 3,
    color = "grey30"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin = margin(t = 3, r = 4, b = 1, l = 4)
  )

# Heat map ONLY (legend removed completely)
pJ_map_only <- ggplot() +
  geom_spatraster(data = heat_factor_2050) +
  scale_fill_manual(
    values = heat_cols,
    guide = "none",
    na.translate = FALSE,
    drop = FALSE
  ) +
  theme_minimal(base_size = 8.2) %+replace%
  theme(
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    axis.title = element_blank(),
    axis.text = element_text(size = 6.5, color = "grey40"),
    panel.grid = element_line(color = "grey93", linewidth = 0.2),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.35),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  )

# Fully custom two-row legend built as text + color boxes
pJ_legend <- ggplot() +
  # Legend title
  annotate(
    "text",
    x = 0.00, y = 1.00,
    label = "2050 Summer Heat Class",
    hjust = 0, vjust = 1,
    size = 5,
    fontface = "bold",
    color = "grey15"
  ) +
  # Row 1, item 1
  annotate("rect", xmin = 0.00, xmax = 0.045, ymin = 0.56, ymax = 0.76,
           fill = "#0072B2", color = "grey30", linewidth = 0.2) +
  annotate("text", x = 0.06, y = 0.66,
           label = "<= 38°C: tolerance",
           hjust = 0, vjust = 0.5, size = 4, color = "grey25") +
  # Row 1, item 2
  annotate("rect", xmin = 0.46, xmax = 0.505, ymin = 0.56, ymax = 0.76,
           fill = "#E69F00", color = "grey30", linewidth = 0.2) +
  annotate("text", x = 0.52, y = 0.66,
           label = "38-40°C: heat stress",
           hjust = 0, vjust = 0.5, size = 4, color = "grey25") +
  # Row 2, item 3
  annotate("rect", xmin = 0.00, xmax = 0.045, ymin = 0.18, ymax = 0.38,
           fill = "#D55E00", color = "grey30", linewidth = 0.2) +
  annotate("text", x = 0.06, y = 0.28,
           label = "> 40°C:  unsuitable",
           hjust = 0, vjust = 0.5, size = 4, color = "grey25") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin = margin(t = 0, r = 4, b = 1, l = 4)
  )

# Caption block
heat_caption <- paste(
  "Physiological thermal limits:",
  "• Blue (<= 38°C): Within tolerance",
  "• Orange (38-40°C): Heat stress / lower activity",
  "• Vermillion (> 40°C): Severe thermal limit",
  sep = "\n"
)

pJ_caption <- ggplot() +
  annotate(
    "text",
    x = 0, y = 1,
    label = heat_caption,
    hjust = 0, vjust = 1,
    size = 4,
    fontface = "bold",
    color = "grey20",
    lineheight = 1.08
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin = margin(t = 1, r = 4, b = 3, l = 4)
  )

# Assemble Panel J
pJ_inner <- cowplot::plot_grid(
  pJ_title,
  pJ_map_only,
  pJ_legend,
  pJ_caption,
  ncol = 1,
  rel_heights = c(0.85, 5.8, 1.2, 1.3)
)

pJ_heat_2050 <- cowplot::ggdraw() +
  cowplot::draw_plot(
    ggplot() +
      theme_void() +
      theme(
        plot.background = element_rect(
          fill = "#FAFAFA",
          color = "grey75",
          linewidth = 0.6
        ),
        plot.margin = margin(0, 0, 0, 0)
      ),
    x = 0, y = 0, width = 1, height = 1
  ) +
  cowplot::draw_plot(
    pJ_inner,
    x = 0.02, y = 0.02,
    width = 0.96, height = 0.96
  )

# =========================================================================
# 5. Final Composite
# =========================================================================
final_figure_4 <- (figure_3x3_clean | pJ_heat_2050) +
  plot_layout(widths = c(3.10, 1.30)) +
  plot_annotation(
    title = "Figure 5. Spatio-temporal trajectories of Parthenium hysterophorus and Zygogramma bicolorata suitability with 2050 summer heat stress",
    subtitle = sprintf(
      "Panels A–I: India-only MaxEnt model and beetle thermal proxy across climate horizons (Parthenium threshold = %.4f, Beetle threshold = %.2f). Panel J: 2050 upper thermal limit constraint.",
      p_thresh,
      bio_thresh
    ),
    theme = theme(
      plot.title = element_text(
        size = 13.5,
        face = "bold",
        color = "grey15",
        hjust = 0,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 9.2,
        color = "grey35",
        hjust = 0,
        margin = margin(b = 6)
      ),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  )

ggsave(
  filename = file.path(out_dir, "Figure4_landscape_full.png"),
  plot = final_figure_4,
  width = 16.5,
  height = 10.5,
  scale = 1.15,
  dpi = 300,
  bg = "white"
)
