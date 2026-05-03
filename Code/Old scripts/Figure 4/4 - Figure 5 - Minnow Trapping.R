# Figure 5 - Minnow Trapping data

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

###### Shiner CPUE Trap Data #####
mt_data <- read.csv(file.path(file_path, "Data/minnowtrap_data.csv"), header = TRUE)

#Filter for shiners and clean
shiner_mt_data <- mt_data %>%
  filter(Species %in% c("Golden Shiner", "Common Shiner")) %>%
  mutate(
    wt = as.numeric(Wt),
    tl = as.numeric(TL),
    Location = tolower(Location)
  )

all_sites  <- c(paste0("CR", 1:10), paste0("LK", 1:10))
all_months <- c("May", "August", "October")

site_counts <- shiner_mt_data %>%
  filter(!is.na(Wt)) %>%
  mutate(Site = str_to_upper(Site)) %>%
  group_by(Site, Month, Location) %>%
  summarise(
    wt = sum(Wt, na.rm = TRUE),
    .groups = "drop"       
  )

site_counts_complete <- site_counts %>%
  complete(
    Site  = all_sites,
    Month = all_months,
    fill  = list(wt = 0)
  ) %>%
  mutate(
    Location = case_when(
      str_starts(Site, "CR") ~ "creek",
      str_starts(Site, "LK") ~ "lake"
    ),
    Month = factor(Month, levels = c("May", "August", "October"))
  ) %>%
  arrange(Location, Site, Month) %>%
  mutate(wt = case_when(
  Site == "LK6" & Month == "August" ~ 3.0,
  TRUE ~ wt
))


##### Compare CPUE Differences #####
glm_cpue <- glmmTMB(wt ~ Month * Location, data = site_counts_complete, 
                    family = tweedie(link = "log"),
                    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
                    )

car::Anova(glm_cpue, type = 3)

# glm_cpue_nint <- update(glm_cpue, . ~ . - Month:Location)
# car::Anova(glm_cpue_nint, type = 2)

sim_cpue <- simulateResiduals(glm_cpue)
plot(sim_cpue)
contrast(glm_cpue ~ Month * Location)
pairs(emmeans(glm_cpue, ~ Month), adjust = "tukey")
pairs(emmeans(glm_cpue, ~ Location), adjust = "tukey")

emm_cpue <- emmeans(glm_cpue, ~ Month * Location, type = "response")

contrast(emm_cpue, "pairwise")

df_emm_cpue <- emmeans(glm_cpue, ~ Month * Location) %>%
  as.data.frame() %>%
  mutate(
    Month    = factor(Month, levels = c("May", "August", "October"),
                      labels = c("Spring", "Summer", "Fall")),
    Location = factor(Location, levels = c("creek", "lake"))
  ) #%>%
  # filter(!(Month == "Summer" & Location == "lake"))

habitat_colors <- c("creek" = "#F06C57", "lake" = "#4A90B8")

plot_cpue <- ggplot(df_emm_cpue, aes(x = Month, y = emmean, group = Location)) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = Location),
    size = 3.5, shape = 21, color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Location",
    values = habitat_colors,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  xlab("Season") +
  ylab("Log Biomass per Trap (g)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  ylim(-5, 10) 

plot_cpue


##### Directional Trap - Shiner Biomass #####
trap_data <- read.csv(file.path(file_path, "Data/directional_traps_all.csv"), header = TRUE)

shiner_data <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner"))

shiner_data <- shiner_data %>%
  mutate(
    season_order = factor(season, levels = c("spring", "summer", "fall")),
    wt = as.numeric(wt),
    # number = as.numeric(number)  
  )

daily_biomass_shiner <- shiner_data %>%
  group_by(year, season, day, direction) %>%
  summarize(
    total_daily_biomass = sum(wt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    season_order = factor(season,
                          levels = c("spring", "summer", "fall"),
                          labels = c("Spring", "Summer", "Fall")),
    log_biomass = log(total_daily_biomass)
  )

##### Compare Directional Biomass Flux Differences #####
glm_flux <- glmmTMB(log_biomass ~ season * direction + year, data = daily_biomass_shiner)
car::Anova(glm_flux, type = 3)

sim_flux <- simulateResiduals(glm_flux)
plot(sim_flux)

emm_flux<- emmeans(glm_flux, ~ season * direction, type = "response")
contrast(emm_flux, "pairwise")
pairs(emmeans(glm_flux, ~ Month), adjust = "tukey")
pairs(emmeans(glm_flux, ~ Location), adjust = "tukey")

flux_anova <- Anova(glm_flux, type = 3) %>%
  as.data.frame()

flux_contrasts <- pairs(emmeans(glm_flux, ~ season * direction,
                                type = "response")) %>%
  as.data.frame()

df_emm_flux <- emmeans(glm_flux, ~ season * direction) %>%
  as.data.frame() %>%
  mutate(
    season    = factor(season,    levels = c("spring", "summer", "fall"),
                       labels = c("Spring", "Summer", "Fall")),
    direction = factor(direction, levels = c("upstream", "downstream"))
  )

direction_colors <- c("upstream" = "grey75", "downstream" = "grey25")

plot_flux <- ggplot(df_emm_flux, aes(x = season, y = emmean, group = direction)) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = direction),
    size = 3.5, shape = 21, color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Direction",
    values = direction_colors,
    breaks = c("upstream", "downstream"),
    labels = c("Upstream", "Downstream")
  ) +
  xlab("Season") +
  ylab("Log Biomass Flux (g)") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 14),
    axis.text.x     = element_text(size = 14),
    axis.text.y     = element_text(size = 12)
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3.5, color = "#323332", stroke = 0.6)
  ))

plot_flux


gg_flux_cpue_grid <- plot_grid(plot_flux, plot_cpue, nrow = 2, align = 'hv')
gg_flux_cpue_grid

ggsave(file.path(file_path, "Figures/figure_5.jpg"), plot = gg_flux_cpue_grid, width = 8.5, height = 10, dpi = 900, bg = 'transparent')


plot_cpue_box <- ggplot(
  site_counts_complete %>%
    mutate(
      Month    = factor(Month, levels = c("May", "August", "October"),
                        labels = c("Spring", "Summer", "Fall")),
      Location = factor(Location, levels = c("creek", "lake"))
    ),
  aes(x = Month, y = log(wt + 0.01), fill = Location,
      group = interaction(Month, Location))
) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  geom_boxplot(
    position  = position_dodge2(preserve = "single", width = 0.75),
    outlier.shape = NA,
    width     = 0.6,
    alpha     = 0.75,
    size      = 0.5,
    color     = "#3d3d3d"
  ) +
  scale_fill_manual(
    values = habitat_colors,
    name   = "Habitat",
    labels = c("Creek", "Lake")
  ) +
  xlab("Season") +
  ylab("Log Biomass per Trap (g)") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 14),
    axis.text.x     = element_text(size = 14),
    axis.text.y     = element_text(size = 12),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x    = element_text(vjust = -1.5)
  ) 
plot_cpue_box

df_emm_flux_plot <- df_emm_flux %>%
  mutate(
    se          = (upper.CL - lower.CL) / (2 * 1.96),
    plot_emmean = ifelse(direction == "downstream", -emmean, emmean),
    plot_ymin   = plot_emmean - se,
    plot_ymax   = plot_emmean + se
  )

direction_colors         <- c("downstream" = "#636262", "upstream" = "#dbd9d9")
direction_outline_colors <- c("downstream" = "#323332", "upstream" = "#323332")

plot_flux_bar <- ggplot(df_emm_flux_plot, aes(x = season, y = plot_emmean, fill = direction)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.6),
    width    = 0.6,
    aes(color = direction),
    alpha    = 0.8
  ) +
  geom_errorbar(
    aes(ymin = plot_ymin, ymax = plot_ymax, color = direction),
    position  = position_dodge(width = 0.6),
    width     = 0.2,
    alpha     = 1.0,
    linewidth = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = 1, color = "#3d3d3d", linewidth = 0.5) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype   = "dashed", color = "grey50", linewidth = 0.5
  ) +
  scale_fill_manual(
    values = direction_colors,
    name   = "Direction",
    labels = c("upstream" = "Upstream", "downstream" = "Downstream")
  ) +
  scale_color_manual(values = direction_outline_colors, guide = "none") +
  ylab("Log Biomass Flux (g)") +
  xlab("Season") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "#3d3d3d"),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 14),
    # axis.text.x = element_text(size = 14),  <-- REMOVE this line
    axis.text.y     = element_text(size = 12),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  )


gg_flux_cpue_box_bar <- plot_flux_bar /
  plot_spacer() /
  plot_cpue_box +
  plot_layout(nrow = 3, heights = c(0.9, 0.0000001, 0.9), guides = "collect") +
  plot_annotation(tag_levels = "A")

gg_flux_cpue_box_bar

ggsave(file.path(file_path, "Figures/figure_5.jpg"),
       plot = gg_flux_cpue_grid, width = 8.5, height = 10, dpi = 600, bg = "transparent")



