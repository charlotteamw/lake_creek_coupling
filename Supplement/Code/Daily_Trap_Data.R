# Supplementary figure for daily trap data

# Author(s): Charlotte Ward & Reilly O'Connor
# Version: 2026-03-30

# Load Pkgs
library(tidyverse)
library(patchwork)
library(grid)
library(car)
library(easystats)
library(cowplot)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(emmeans)

file_path <- getwd()

source(file.path(file_path, "/Code/0 - Functions.R"))


##### Directional Trap - Shiner Biomass (Weir & Fyke, 2023 & 2024) #####
trap_data <- read.csv(file.path(file_path, "Data/directional_traps_all.csv"), header = TRUE)

shiner_data <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner")) %>%
  mutate(
    wt        = as.numeric(wt),
    month_low = tolower(month),
    direction = factor(direction, levels = c("upstream", "downstream"))
  ) %>%
  # Remove April 24, 2024 from all trap types
  filter(!(year == 2024 & month_low == "april" & day == 24))

direction_colors <- c("upstream" = "grey75", "downstream" = "grey25")

# --- Build daily biomass summary by trap type and year ---
make_daily_biomass <- function(df, trap_type, yr) {
  df %>%
    filter(year == yr, type == trap_type) %>%
    mutate(
      # Reclassify late-September days within the "october" sampling period
      true_month = case_when(
        month_low == "october" & day >= 25 ~ "september",
        TRUE ~ month_low
      ),
      date = as.Date(paste(yr, true_month, day), format = "%Y %B %d"),
      season = case_when(
        month_low %in% c("april", "may") ~ "Spring",
        month_low == "august"            ~ "Summer",
        month_low == "october"           ~ "Fall"
      ),
      season = factor(season, levels = c("Spring", "Summer", "Fall"))
    ) %>%
    group_by(season, date, direction) %>%
    summarize(
      total_daily_biomass = sum(wt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(direction = factor(direction, levels = c("upstream", "downstream")))
}

# --- Build the faceted point + line plot ---
make_daily_plot <- function(df, trap_label, yr) {
  ggplot(df, aes(x = date, y = total_daily_biomass,
                 color = direction, fill = direction, group = direction)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 3, shape = 21, color = "#323332", stroke = 0.5) +
    scale_fill_manual(
      name   = "Direction",
      values = direction_colors,
      labels = c("Upstream", "Downstream")
    ) +
    scale_color_manual(
      name   = "Direction",
      values = direction_colors,
      labels = c("Upstream", "Downstream")
    ) +
    scale_x_date(
      date_labels = "%b %d",
      date_breaks = "2 days"
    ) +
    facet_wrap(~ season, nrow = 1, scales = "free_x") +
    labs(
      title = paste0("Directional Trap \u2013 ", trap_label, " (", yr, ")"),
      x     = "Date",
      y     = "Total Daily Biomass (g)"
    ) +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black"),
      strip.text       = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "black"),
      legend.position  = "right",
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 11),
      axis.title       = element_text(size = 13),
      axis.text.x      = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 11),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, size = 3,
                                               color = "#323332", stroke = 0.5)),
      color = "none"
    )
}

# --- Build weir plots ---
daily_weir_2023 <- make_daily_biomass(shiner_data, "weir", 2023)
daily_weir_2024 <- make_daily_biomass(shiner_data, "weir", 2024)

plot_weir_2023  <- make_daily_plot(daily_weir_2023, "Weir", 2023)
plot_weir_2024  <- make_daily_plot(daily_weir_2024, "Weir", 2024)

# --- Build fyke plots ---
daily_fyke_2023 <- make_daily_biomass(shiner_data, "fyke", 2023)
daily_fyke_2024 <- make_daily_biomass(shiner_data, "fyke", 2024)

plot_fyke_2023  <- make_daily_plot(daily_fyke_2023, "Fyke", 2023)
plot_fyke_2024  <- make_daily_plot(daily_fyke_2024, "Fyke", 2024)

# --- Combine and save: Weir (2023 upper, 2024 lower) ---
gg_weir_grid <- plot_grid(plot_weir_2023, plot_weir_2024, nrow = 2, align = "hv")
gg_weir_grid

ggsave(file.path(file_path, "Figures/figure_5_weir_2023_2024.jpg"),
       plot = gg_weir_grid, width = 13, height = 8, dpi = 600, bg = "transparent")

# --- Combine and save: Fyke (2023 upper, 2024 lower) ---
gg_fyke_grid <- plot_grid(plot_fyke_2023, plot_fyke_2024, nrow = 2, align = "hv")
gg_fyke_grid

ggsave(file.path(file_path, "Figures/figure_5_fyke_2023_2024.jpg"),
       plot = gg_fyke_grid, width = 13, height = 8, dpi = 600, bg = "transparent")
