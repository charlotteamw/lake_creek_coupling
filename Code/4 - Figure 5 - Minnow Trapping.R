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
    size = 6, shape = 21, color = "#323332", stroke = 0.6,
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
  group_by(season, day, direction) %>%
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
glm_flux <- glmmTMB(log_biomass ~ season * direction, data = daily_biomass_shiner)
car::Anova(glm_flux, type = 3)

sim_flux <- simulateResiduals(glm_flux)
plot(sim_flux)

emm_flux<- emmeans(glm_flux, ~ season * direction, type = "response")
contrast(emm_flux, "pairwise")
pairs(emmeans(glm_flux, ~ Month), adjust = "tukey")
pairs(emmeans(glm_flux, ~ Location), adjust = "tukey")

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
    size = 6, shape = 21, color = "#323332", stroke = 0.6,
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
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 6, color = "#323332", stroke = 0.6)
  ))

plot_flux


gg_flux_cpue_grid <- plot_grid(plot_flux, plot_cpue, nrow = 2, align = 'hv')
gg_flux_cpue_grid

ggsave(file.path(file_path, "Figures/figure_5.jpg"), plot = gg_flux_cpue_grid, width = 8.5, height = 10, dpi = 600, bg = 'transparent')


