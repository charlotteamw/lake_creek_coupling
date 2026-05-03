# Figure 5 - Directional Trap Flux + Daily Weir/Fyke 2023 & 2024
# Author(s): Charlotte Ward & Reilly O'Connor
# Version: 2026-04-09

# ── Packages ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(cowplot)
library(car)
library(easystats)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(emmeans)

file_path <- getwd()
source(file.path(file_path, "/Code/0 - Functions.R"))

direction_colors <- c("upstream" = "grey75", "downstream" = "grey25")

# ─────────────────────────────────────────────────────────────────────────────
# PANEL A: Directional Trap Flux (emmeans)
# ─────────────────────────────────────────────────────────────────────────────
trap_data <- read.csv(file.path(file_path, "Data/directional_traps_all.csv"),
                      header = TRUE)

shiner_data_flux <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner")) %>%
  mutate(
    wt           = as.numeric(wt),
    season_order = factor(season, levels = c("spring", "summer", "fall"))
  )

daily_biomass_shiner <- shiner_data_flux %>%
  group_by(season, day, direction) %>%
  summarize(total_daily_biomass = sum(wt, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    season_order = factor(season,
                          levels = c("spring", "summer", "fall"),
                          labels = c("Spring", "Summer", "Fall")),
    log_biomass  = log(total_daily_biomass)
  )

glm_flux <- glmmTMB(log_biomass ~ season * direction,
                    data = daily_biomass_shiner)
car::Anova(glm_flux, type = 3)

sim_flux <- simulateResiduals(glm_flux)
plot(sim_flux)

emm_flux <- emmeans(glm_flux, ~ season * direction, type = "response")
contrast(emm_flux, "pairwise")

df_emm_flux <- emmeans(glm_flux, ~ season * direction) %>%
  as.data.frame() %>%
  mutate(
    season    = factor(season,
                       levels = c("spring", "summer", "fall"),
                       labels = c("Spring", "Summer", "Fall")),
    direction = factor(direction, levels = c("upstream", "downstream"))
  )

# ── Plot A — add title, remove legend ────────────────────────────────────────
plot_flux <- ggplot(df_emm_flux,
                    aes(x = season, y = emmean, group = direction)) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = direction),
    size = 3, shape = 21, color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Direction",
    values = direction_colors,
    breaks = c("upstream", "downstream"),
    labels = c("Upstream", "Downstream")
  ) +
  geom_vline(xintercept = c(1.5, 2.5),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(title = "Seasonal Summary", x = NULL, y = "Log Daily Biomass (g)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.line       = element_line(color = "black"),
    legend.position = "none",                # legend removed here
    axis.title      = element_text(size = 11),
    axis.text.x     = element_text(size = 11),
    axis.text.y     = element_text(size = 10),
    plot.title      = element_text(size = 11),
  )

plot_flux
# ── Updated make_daily_plot — accepts custom title, legend removed ────────────
make_daily_plot <- function(df, plot_title) {
  
  all_dates       <- sort(unique(df$date))
  date_labels     <- format(all_dates, "%b %d")
  
  df <- df %>%
    mutate(date_f = factor(format(date, "%b %d"), levels = date_labels))
  
  spring_last_idx <- max(which(all_dates %in% df$date[df$season == "Spring"]))
  summer_last_idx <- max(which(all_dates %in% df$date[df$season == "Summer"]))
  vline_positions <- c(spring_last_idx + 0.5, summer_last_idx + 0.5)
  
  ggplot(df, aes(x = date_f, y = total_daily_biomass,
                 color = direction, fill = direction,
                 group = interaction(direction, season))) +
    geom_vline(xintercept = vline_positions,
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.5, shape = 21, color = "#323332", stroke = 0.6) +
    scale_fill_manual(
      name   = "Direction",
      values = direction_colors,
      breaks = c("upstream", "downstream"),
      labels = c("Upstream", "Downstream")
    ) +
    scale_color_manual(
      name   = "Direction",
      values = direction_colors,
      breaks = c("upstream", "downstream"),
      labels = c("Upstream", "Downstream")
    ) +
    scale_x_discrete(labels = date_labels) +
    labs(title = plot_title, x = NULL, y = "Daily Biomass (g)") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.line       = element_line(color = "black"),
      legend.position = "none",                # legend removed here
      axis.title      = element_text(size = 11),
      axis.text.x     = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 10),
      plot.title      = element_text(size = 11)
    )
}

# ── Build all four daily plots with custom titles ─────────────────────────────
plot_weir_2023 <- make_daily_plot(daily_weir_2023, "Culvert 2023")
plot_fyke_2023 <- make_daily_plot(daily_fyke_2023, "Outflow 2023")
plot_weir_2024 <- make_daily_plot(daily_weir_2024, "Culvert 2024")
plot_fyke_2024 <- make_daily_plot(daily_fyke_2024, "Outflow 2024")

# ── Extract shared legend from a plot that still has one ─────────────────────
legend_plot <- ggplot(df_emm_flux,
                      aes(x = season, y = emmean, group = direction)) +
  geom_point(aes(fill = direction), size = 3.5, shape = 21,
             color = "#323332", stroke = 0.6) +
  scale_fill_manual(
    name   = "Direction",
    values = direction_colors,
    breaks = c("upstream", "downstream"),
    labels = c("Upstream", "Downstream")
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 10)
  ) +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3.5, color = "#323332", stroke = 0.6)
  ))

shared_legend <- get_legend(legend_plot)

# ── Combine panels (no individual legends) + shared legend on right ───────────


plots_no_legend <- plot_grid(
  plot_flux,
  plot_fyke_2023,
  plot_weir_2023,
  plot_fyke_2024,
  plot_weir_2024,
  nrow       = 5,
  align      = "hv",
  greedy     = FALSE
)

gg_combined <- plot_grid(
  plots_no_legend,
  shared_legend,
  ncol        = 2,
  rel_widths  = c(1, 0.18)   # narrow right column for legend
)

gg_combined

ggsave(
  file.path(file_path, "Figures/figure_5_flux_daily_all.jpg"),
  plot   = gg_combined,
  width  = 8,
  height = 9,
  dpi    = 900,
  bg     = "transparent"
)
