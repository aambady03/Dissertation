library(ggplot2)
library(patchwork)
library(terra)
library(tidyterra)
library(grid)
library(cowplot)

# =========================================================================
# 1. Base Theme — larger text, no subtitles
# =========================================================================
theme_landscape_map <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title        = element_text(size = rel(1.0), face = "bold",
                                       color = "grey15", hjust = 0.5,
                                       margin = margin(b = 3)),
      plot.subtitle     = element_blank(),   # removed from all panels
      axis.title        = element_blank(),
      axis.text         = element_text(size = rel(0.75), color = "grey40"),
      panel.grid        = element_line(color = "grey93", linewidth = 0.2),
      panel.border      = element_rect(color = "grey80", fill = NA,
                                       linewidth = 0.35),
      legend.title      = element_text(size = rel(0.90), face = "bold",
                                       color = "grey20"),
      legend.text       = element_text(size = rel(0.82), color = "grey30"),
      legend.position   = "right",
      legend.key.height = unit(0.55, "cm"),
      legend.key.width  = unit(0.38, "cm"),
      legend.margin     = margin(l = 3, r = 3),
      plot.margin       = margin(t = 3, r = 3, b = 3, l = 3)
    )
}

# =========================================================================
# 2. Plotting Functions
#    - Suitability: legend REMOVED (shared legend added via patchwork)
#    - Overlap:     legend REMOVED (shared legend added via patchwork)
# =========================================================================

# --- Suitability map (no legend — shared per row) ---
plot_suitability_noleg <- function(r, title_str) {
  ggplot() +
    geom_spatraster(data = r) +
    scale_fill_viridis_c(
      option    = "C",
      limits    = c(0, 1),
      name      = "Suitability\n(cloglog)",
      na.value  = "transparent",
      guide     = "none"          # legend suppressed here
    ) +
    labs(title = title_str) +
    theme_landscape_map()
}

# --- Suitability map WITH legend (rightmost column only) ---
plot_suitability_leg <- function(r, title_str) {
  ggplot() +
    geom_spatraster(data = r) +
    scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      name   = "Suitability\n(cloglog)",
      na.value = "transparent"
    ) +
    labs(title = title_str) +
    guides(fill = guide_colorbar(
      barheight       = unit(4.0, "cm"),
      barwidth        = unit(0.45, "cm"),
      ticks.linewidth = 0.4,
      frame.colour    = "grey50",
      frame.linewidth = 0.4
    )) +
    theme_landscape_map()
}

# --- Overlap map (no legend) ---
plot_overlap_noleg <- function(fact_r, title_str) {
  overlap_cols <- c(
    "Neither suitable" = "#D9D9D9",
    "Parthenium only"  = "#E69F00",
    "Beetle only"      = "#0072B2",
    "Both suitable"    = "#009E73"
  )
  ggplot() +
    geom_spatraster(data = fact_r) +
    scale_fill_manual(
      values       = overlap_cols,
      name         = "Categorical\nOverlap",
      na.translate = FALSE,
      drop         = FALSE,
      guide        = "none"       # legend suppressed here
    ) +
    labs(title = title_str) +
    theme_landscape_map()
}

# --- Overlap map WITH legend (rightmost column only) ---
plot_overlap_leg <- function(fact_r, title_str) {
  overlap_cols <- c(
    "Neither suitable" = "#D9D9D9",
    "Parthenium only"  = "#E69F00",
    "Beetle only"      = "#0072B2",
    "Both suitable"    = "#009E73"
  )
  ggplot() +
    geom_spatraster(data = fact_r) +
    scale_fill_manual(
      values       = overlap_cols,
      name         = "Categorical\nOverlap",
      na.translate = FALSE,
      drop         = FALSE
    ) +
    labs(title = title_str) +
    guides(fill = guide_legend(
      keyheight    = unit(0.55, "cm"),
      keywidth     = unit(0.55, "cm"),
      override.aes = list(color = "grey40", linewidth = 0.3)
    )) +
    theme_landscape_map()
}

# =========================================================================
# 3. Column header strips
#    These sit above each column as plain text panels
# =========================================================================
make_col_header <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label    = label,
      hjust    = 0.5, vjust = 0.5,
      size     = 5.2,
      fontface = "bold",
      color    = "grey10"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(t = 4, b = 2, l = 2, r = 2))
}

hdr_current <- make_col_header("Current")
hdr_2030    <- make_col_header("2030 (SSP2-4.5)")
hdr_2050    <- make_col_header("2050 (SSP2-4.5)")

# Empty top-left corner to align headers with the 3-column grid
hdr_blank   <- ggplot() + theme_void()

# =========================================================================
# 4. Row labels (left side)
# =========================================================================
make_row_label <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label    = label,
      hjust    = 0.5, vjust = 0.5,
      size     = 4.2,
      fontface = "bold",
      color    = "grey15",
      angle    = 90
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(t = 2, b = 2, l = 2, r = 2))
}

row_lbl_parth  <- make_row_label("Parthenium suitability")
row_lbl_beetle <- make_row_label("Beetle thermal proxy")
row_lbl_over   <- make_row_label("Niche overlap")

# =========================================================================
# 5. Build the 9 map panels
#    Columns 1 & 2 = no legend | Column 3 = legend (shared per row)
# =========================================================================

# --- Row A: Parthenium suitability ---
pA <- plot_suitability_noleg(pair_now$parthenium, "A.")
pB <- plot_suitability_noleg(pair_30$parthenium,  "B.")
pC <- plot_suitability_leg( pair_50$parthenium,  "C.")

# --- Row B: Beetle thermal proxy ---
pD <- plot_suitability_noleg(pair_now$beetle, "D.")
pE <- plot_suitability_noleg(pair_30$beetle,  "E.")
pF <- plot_suitability_leg( pair_50$beetle,  "F.")

# --- Row C: Niche overlap ---
pG <- plot_overlap_noleg(pair_now$overlap_factor, "G.")
pH <- plot_overlap_noleg(pair_30$overlap_factor,  "H.")
pI <- plot_overlap_leg(  pair_50$overlap_factor,  "I.")

# =========================================================================
# 6. Assemble 3x3 grid with column headers and row labels
#
#   Layout (using patchwork design string):
#   [blank] [hdr_cur] [hdr_30] [hdr_50]
#   [rowA ] [pA     ] [pB    ] [pC    ]
#   [rowB ] [pD     ] [pE    ] [pF    ]
#   [rowC ] [pG     ] [pH    ] [pI    ]
# =========================================================================
grid_3x3 <- (
  # Row 0: headers
  (hdr_blank | hdr_current | hdr_2030 | hdr_2050) /
    
    # Row 1: Parthenium
    (row_lbl_parth | pA | pB | pC) /
    
    # Row 2: Beetle
    (row_lbl_beetle | pD | pE | pF) /
    
    # Row 3: Overlap
    (row_lbl_over | pG | pH | pI)
) +
  plot_layout(
    heights = c(0.08, 1, 1, 1),   # header row is compact
    widths  = c(0.07, 1, 1, 1)    # row label column is compact
  )

# =========================================================================
# 7. Panel J — 2050 Summer Heat Stress
#    Keep only the thermal limits caption, remove all other captions
# =========================================================================
heat_cols <- c(
  "<= 38°C: within tolerance" = "#0072B2",
  "38–40°C: heat stress"      = "#E69F00",
  "> 40°C: likely unsuitable" = "#D55E00"
)

# Map only — no legend
pJ_map_only <- ggplot() +
  geom_spatraster(data = heat_factor_2050) +
  scale_fill_manual(
    values       = heat_cols,
    guide        = "none",
    na.translate = FALSE,
    drop         = FALSE
  ) +
  labs(title = "J. 2050 Summer Heat Stress") +
  theme_landscape_map() +
  theme(
    plot.title      = element_text(
      size = 11, face = "bold",
      color = "#8B0000", hjust = 0
    ),
    plot.background = element_rect(fill = "#FAFAFA", color = NA)
  )

# Custom legend panel
pJ_legend <- ggplot() +
  annotate("text",
           x = 0.00, y = 0.96,
           label    = "Thermal limit class",
           hjust    = 0, vjust = 1,
           size     = 4.2, fontface = "bold", color = "grey15") +
  # Item 1
  annotate("rect",
           xmin = 0.00, xmax = 0.07, ymin = 0.70, ymax = 0.84,
           fill = "#0072B2", color = "grey30", linewidth = 0.25) +
  annotate("text",
           x = 0.11, y = 0.77,
           label = "\u2264 38\u00b0C: within tolerance",
           hjust = 0, vjust = 0.5, size = 3.8, color = "grey20") +
  # Item 2
  annotate("rect",
           xmin = 0.00, xmax = 0.07, ymin = 0.50, ymax = 0.64,
           fill = "#E69F00", color = "grey30", linewidth = 0.25) +
  annotate("text",
           x = 0.11, y = 0.57,
           label = "38\u201340\u00b0C: heat stress",
           hjust = 0, vjust = 0.5, size = 3.8, color = "grey20") +
  # Item 3
  annotate("rect",
           xmin = 0.00, xmax = 0.07, ymin = 0.30, ymax = 0.44,
           fill = "#D55E00", color = "grey30", linewidth = 0.25) +
  annotate("text",
           x = 0.11, y = 0.37,
           label = "> 40\u00b0C: likely unsuitable",
           hjust = 0, vjust = 0.5, size = 3.8, color = "grey20") +
  # Thermal limits caption (KEPT as requested)
  annotate("text",
           x = 0.00, y = 0.14,
           label = "Physiological upper limits\nfor Z. bicolorata activity",
           hjust = 0, vjust = 1,
           size = 3.4, color = "grey40", lineheight = 1.1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin     = margin(t = 4, r = 6, b = 4, l = 4)
  )

# Assemble Panel J
pJ_inner <- cowplot::plot_grid(
  pJ_map_only,
  pJ_legend,
  ncol        = 1,
  rel_heights = c(4.5, 1.8)
)

pJ_heat_2050 <- cowplot::ggdraw() +
  cowplot::draw_plot(
    ggplot() +
      theme_void() +
      theme(
        plot.background = element_rect(
          fill      = "#FAFAFA",
          color     = "grey75",
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
# 8. Final Composite
# =========================================================================
final_figure_4 <- (grid_3x3 | pJ_heat_2050) +
  plot_layout(widths = c(3.10, 1.10)) +
  plot_annotation(
    title = paste0(
      "Figure 5. Spatio-temporal trajectories of Parthenium hysterophorus ",
      "and Zygogramma bicolorata suitability with 2050 summer heat stress"
    ),
    theme = theme(
      plot.title = element_text(
        size = 13.5, face = "bold", color = "grey15",
        hjust = 0, margin = margin(b = 3)
      ),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  )

ggsave(
  filename = file.path(out_dir, "Figure5_landscape_full.png"),
  plot     = final_figure_4,
  width    = 16.5,
  height   = 10.5,
  scale    = 1.15,
  dpi      = 300,
  bg       = "white"
)

cat("✅ Saved: Figure4_landscape_full.png\n")