############################################################
# SUPPLEMENTARY BIPLOTS
# Plot 1: Shiner individuals + baseline means ± SD
# Plot 2: All consumer species means ± SD + baseline means ± SD
############################################################

# ==============================
# Load packages
# ==============================
library(tidyverse)
library(ggrepel)

# ==============================
# Helpers & shared aesthetics
# ==============================
lake_creek_colors <- c("Creek" = "#F06C57", "Lake" = "#4A90B8")

custom_theme <- theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    axis.title   = element_text(size = 13),
    axis.line    = element_line(color = "black"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )

# ==============================
# Load data
# ==============================
isotope_data <- read.csv(
  "~/Documents/lake_creek_coupling/Data/iso_metadata.csv",
  header = TRUE
)

# ==============================
# Shared derived objects
# ==============================

# Baseline means ± SD (mayfly & mussel), title-cased location
baseline_means_sd <- isotope_data %>%
  filter(organism %in% c("mayfly", "mussel")) %>%
  mutate(
    organism = stringr::str_to_title(organism),
    location = stringr::str_to_title(location)
  ) %>%
  group_by(location, organism) %>%
  summarise(
    d13C_mean = mean(d13C, na.rm = TRUE),
    d13C_sd   = sd(d13C,   na.rm = TRUE),
    d15N_mean = mean(d15N, na.rm = TRUE),
    d15N_sd   = sd(d15N,   na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(label = organism)

# Shiner individual points (liver only), title-cased location
shiner_plot_df <- isotope_data %>%
  filter(
    organism %in% c("golden shiner", "common shiner"),
    tissue == "liver"
  ) %>%
  mutate(location = stringr::str_to_title(location))

# Consumer species means ± SD, title-cased location
abbr <- c(
  "smallmouth bass"      = "SMB",
  "golden shiner"        = "GS",
  "creek chub"           = "CC",
  "common shiner"        = "CS",
  "pumpkinseed"          = "PS",
  "yellow perch"         = "YP",
  "brook trout"          = "BT",
  "pearl dace"           = "PD",
  "blacknose shiner"     = "BNS",
  "northern redbelly dace" = "NRD",
  "bluntnose minnow"     = "BNM"
)

species_means_sd <- isotope_data %>%
  filter(organism %in% names(abbr)) %>%
  mutate(location = stringr::str_to_title(location)) %>%
  group_by(location, organism) %>%
  summarise(
    d13C_mean = mean(d13C, na.rm = TRUE),
    d13C_sd   = sd(d13C,   na.rm = TRUE),
    d15N_mean = mean(d15N, na.rm = TRUE),
    d15N_sd   = sd(d15N,   na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(label_abbr = abbr[organism])

# =============================================================================
# PLOT 1: Shiner individuals + baseline means ± SD
# =============================================================================
plot_biplot_shiners <- ggplot() +
  # Shiner individual points
  geom_point(
    data = shiner_plot_df,
    aes(x = d13C, y = d15N, color = location),
    alpha = 0.6, size = 2
  ) +
  # Baseline vertical error bars (± SD in δ15N)
  geom_errorbar(
    data = baseline_means_sd,
    aes(
      x     = d13C_mean,
      ymin  = d15N_mean - d15N_sd,
      ymax  = d15N_mean + d15N_sd,
      color = location
    ),
    width = 0
  ) +
  # Baseline horizontal error bars (± SD in δ13C)
  geom_segment(
    data = baseline_means_sd,
    aes(
      x     = d13C_mean - d13C_sd,
      xend  = d13C_mean + d13C_sd,
      y     = d15N_mean,
      yend  = d15N_mean,
      color = location
    ),
    linewidth = 0.6
  ) +
  # Baseline mean points (shape by organism)
  geom_point(
    data = baseline_means_sd,
    aes(x = d13C_mean, y = d15N_mean, color = location, shape = organism),
    size = 2.2, stroke = 1.1, fill = "white"
  ) +
  scale_color_manual(values = lake_creek_colors, name = "Location") +
  scale_shape_manual(values = c("Mayfly" = 21, "Mussel" = 24), name = "Baseline") +
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030")
  ) +
  theme_bw() +
  custom_theme

print(plot_biplot_shiners)
ggsave("supp_biplot_shiners.png", plot = plot_biplot_shiners,
       width = 6, height = 5, dpi = 600)

# =============================================================================
# PLOT 2: All consumer species means ± SD + baseline means ± SD
# =============================================================================
plot_biplot_allspecies <- ggplot() +
  # Baseline error bars (lighter, behind consumers)
  geom_errorbar(
    data = baseline_means_sd,
    aes(x = d13C_mean, ymin = d15N_mean - d15N_sd, ymax = d15N_mean + d15N_sd,
        color = location),
    width = 0, linewidth = 0.4, alpha = 0.6
  ) +
  geom_segment(
    data = baseline_means_sd,
    aes(x = d13C_mean - d13C_sd, xend = d13C_mean + d13C_sd,
        y = d15N_mean, yend = d15N_mean, color = location),
    linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = baseline_means_sd,
    aes(x = d13C_mean, y = d15N_mean, color = location),
    shape = 21, size = 2.4, stroke = 0.9, fill = "white"
  ) +
  geom_label_repel(
    data = baseline_means_sd,
    aes(x = d13C_mean, y = d15N_mean, label = label, color = location),
    size = 3.0, label.size = 0, fill = scales::alpha("white", 0.7),
    box.padding = 0.3, point.padding = 0.2, seed = 123, show.legend = FALSE
  ) +
  # Consumer error bars (emphasized, in front)
  geom_errorbar(
    data = species_means_sd,
    aes(x = d13C_mean, ymin = d15N_mean - d15N_sd, ymax = d15N_mean + d15N_sd,
        color = location),
    width = 0, linewidth = 0.6, alpha = 0.75
  ) +
  geom_segment(
    data = species_means_sd,
    aes(x = d13C_mean - d13C_sd, xend = d13C_mean + d13C_sd,
        y = d15N_mean, yend = d15N_mean, color = location),
    linewidth = 0.6, alpha = 0.75
  ) +
  geom_point(
    data = species_means_sd,
    aes(x = d13C_mean, y = d15N_mean, color = location),
    size = 3.2, stroke = 0.7
  ) +
  geom_label_repel(
    data = species_means_sd,
    aes(x = d13C_mean, y = d15N_mean, label = label_abbr, color = location),
    size = 3.2, label.size = 0, segment.size = 0.3,
    fill = scales::alpha("white", 0.65),
    box.padding = 0.4, point.padding = 0.35, max.overlaps = 200,
    seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = lake_creek_colors, name = "Location") +
  labs(
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030")
  ) +
  theme_bw() +
  theme(
    axis.line    = element_line(colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  )

print(plot_biplot_allspecies)
ggsave("supp_biplot_allspecies.png", plot = plot_biplot_allspecies,
       width = 7, height = 5.5, dpi = 600)
