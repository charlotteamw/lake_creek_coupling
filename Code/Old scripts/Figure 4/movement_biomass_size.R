# =============================================================================
# Figure 5 — Flux, CPUE, and Fish Size
# Author(s): Charlotte Ward & Reilly O'Connor
# Updated: 2026-04-14
# Description:
#   Panel 1 = Directional trap flux (emmeans)
#   Panel 2 = Shiner CPUE trap data (emmeans)
#   Panel 3 = Shiner total length (GLM + emmeans)
#   Panel 4 = Smallmouth bass total length (GLM + emmeans; additive model)
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(cowplot)
library(car)
library(easystats)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(emmeans)
library(ggpubr)
library(readr)
library(stringr)
library(tidyr)

# ── Setup ─────────────────────────────────────────────────────────────────────
file_path <- getwd()
source(file.path(file_path, "Code/0 - Functions.R"))

options(contrasts = c("contr.sum", "contr.poly"))

ensure_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

direction_colors <- c("upstream" = "grey75", "downstream" = "grey25")
habitat_colors   <- c("creek" = "#F06C57", "lake" = "#4A90B8")

base_theme <- theme_bw() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    legend.position  = "none",
    axis.title       = element_text(size = 12),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 10)
  )

# =============================================================================
# PANEL 1 — Directional Trap Flux
# =============================================================================

trap_data <- read.csv(
  file.path(file_path, "Data/directional_traps_all.csv"),
  header = TRUE
)

shiner_data_flux <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner")) %>%
  mutate(
    wt = ensure_numeric(wt),
    season_order = factor(season, levels = c("spring", "summer", "fall"))
  )

daily_biomass_shiner <- shiner_data_flux %>%
  group_by(season, day, direction) %>%
  summarise(
    total_daily_biomass = sum(wt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    season_order = factor(
      season,
      levels = c("spring", "summer", "fall"),
      labels = c("Spring", "Summer", "Fall")
    ),
    log_biomass = log(total_daily_biomass)
  )

glm_flux <- glmmTMB(
  log_biomass ~ season * direction,
  data = daily_biomass_shiner
)

car::Anova(glm_flux, type = 3)


flux_anova <- Anova(glm_flux, type = 3) %>%
  as.data.frame()


sim_flux <- simulateResiduals(glm_flux)
plot(sim_flux)

emm_flux <- emmeans(glm_flux, ~ season * direction, type = "response")
contrast(emm_flux, "pairwise")

df_emm_flux <- emmeans(glm_flux, ~ season * direction) %>%
  as.data.frame() %>%
  mutate(
    season = factor(
      season,
      levels = c("spring", "summer", "fall"),
      labels = c("Spring", "Summer", "Fall")
    ),
    direction = factor(direction, levels = c("upstream", "downstream"))
  )

plot_flux <- ggplot(
  df_emm_flux,
  aes(x = season, y = emmean, group = direction)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.6,
    color = "#323332",
    alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = direction),
    size = 3,
    shape = 21,
    color = "#323332",
    stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Direction",
    values = direction_colors,
    breaks = c("upstream", "downstream"),
    labels = c("Upstream", "Downstream")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = NULL,
    y = "Log Daily Biomass (g)"
  ) +
  coord_cartesian(clip = "off") +
  base_theme

# =============================================================================
# PANEL 2 — Shiner CPUE Trap Data
# =============================================================================

mt_data <- read.csv(
  file.path(file_path, "Data/minnowtrap_data.csv"),
  header = TRUE
)

shiner_mt_data <- mt_data %>%
  filter(Species %in% c("Golden Shiner", "Common Shiner")) %>%
  mutate(
    wt = ensure_numeric(Wt),
    tl = ensure_numeric(TL),
    Location = tolower(Location)
  )

all_sites  <- c(paste0("CR", 1:10), paste0("LK", 1:10))
all_months <- c("May", "August", "October")

site_counts <- shiner_mt_data %>%
  filter(!is.na(wt)) %>%
  mutate(Site = str_to_upper(Site)) %>%
  group_by(Site, Month, Location) %>%
  summarise(
    wt = sum(wt, na.rm = TRUE),
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
  mutate(
    wt = case_when(
      Site == "LK6" & Month == "August" ~ 3.0,
      TRUE ~ wt
    )
  )

glm_cpue <- glmmTMB(
  wt ~ Month * Location,
  data = site_counts_complete,
  family = tweedie(link = "log"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000)
  )
)

car::Anova(glm_cpue, type = 3)

cpue_anova <- Anova(glm_cpue, type = 3) %>%
  as.data.frame()

sim_cpue <- simulateResiduals(glm_cpue)
plot(sim_cpue)

emm_cpue <- emmeans(glm_cpue, ~ Month * Location, type = "response")
contrast(emm_cpue, "pairwise")

df_emm_cpue <- emmeans(glm_cpue, ~ Month * Location) %>%
  as.data.frame() %>%
  mutate(
    Month = factor(
      Month,
      levels = c("May", "August", "October"),
      labels = c("Spring", "Summer", "Fall")
    ),
    Location = factor(Location, levels = c("creek", "lake"))
  )

plot_cpue <- ggplot(
  df_emm_cpue,
  aes(x = Month, y = emmean, group = Location)
) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 0.6,
    color = "#323332",
    alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = Location),
    size = 3,
    shape = 21,
    color = "#323332",
    stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Location",
    values = habitat_colors,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = NULL,
    y = "Log Biomass per Trap (g)"
  ) +
  coord_cartesian(clip = "off") +
  base_theme

# =============================================================================
# PANEL 3 — Shiner Size (GLM + EMMs)
# =============================================================================

fish_data <- read.csv(
  file.path(file_path, "Data/fish_catches_all.csv"),
  header = TRUE
)

df_shiner_size <- fish_data %>%
  mutate(
    tl       = ensure_numeric(tl),
    year     = ensure_numeric(year),
    location = str_to_title(location),
    season   = str_to_title(season),
    species  = str_to_lower(species)
  ) %>%
  filter(
    species %in% c("common shiner", "golden shiner"),
    !is.na(tl),
    tl > 0,
    location %in% c("Creek", "Lake"),
    season %in% c("Spring", "Summer", "Fall")
  ) %>%
  mutate(
    season = factor(season, levels = c("Spring", "Summer", "Fall")),
    location = factor(location, levels = c("Creek", "Lake"))
  )

glm_size_shiner <- glm(
  tl ~ season * location,
  data = df_shiner_size,
  family = gaussian()
)

car::Anova(glm_size_shiner, type = 3)


shiner_size_anova <- Anova(glm_size_shiner, type = 3) %>%
  as.data.frame()

sim_size_shiner <- simulateResiduals(glm_size_shiner)
plot(sim_size_shiner)

shiner_size_contrasts <- as.data.frame(
  pairs(emmeans(glm_size_shiner, ~ season * location), adjust = "tukey")
)

par(mfrow = c(2, 2))
plot(glm_size_shiner)


df_emm_size_shiner <- emmeans(glm_size_shiner, ~ season * location) %>%
  as.data.frame() %>%
  mutate(
    season = factor(as.character(season), levels = c("Spring", "Summer", "Fall")),
    location = factor(tolower(as.character(location)), levels = c("creek", "lake"))
  )

plot_size_shiner <- ggplot(
  df_emm_size_shiner,
  aes(x = season, y = emmean, group = location)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.6,
    color = "#323332",
    alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = location),
    size = 3,
    shape = 21,
    color = "#323332",
    stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Location",
    values = habitat_colors,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = NULL,
    y = "Shiner total length (mm)"
  ) +
  coord_cartesian(clip = "off") +
  base_theme


# =============================================================================
# PANEL 4 — Smallmouth Bass Size (Gamma GLM + back-transformed EMMs)
# =============================================================================

df_smb_size <- fish_data %>%
  mutate(
    tl       = ensure_numeric(tl),
    year     = ensure_numeric(year),
    location = str_to_title(location),
    season   = str_to_title(season),
    species  = str_to_lower(species)
  ) %>%
  filter(
    species == "smallmouth bass",
    !is.na(tl),
    tl > 0,
    location %in% c("Creek", "Lake"),
    season %in% c("Spring", "Summer", "Fall")
  ) %>%
  mutate(
    season = factor(season, levels = c("Spring", "Summer", "Fall")),
    location = factor(location, levels = c("Creek", "Lake"))
  )

glm_size_smb_gamma <- glm(
  tl ~ season + location,
  data = df_smb_size,
  family = Gamma(link = "log")
)


car::Anova(glm_size_smb_gamma, type = 2)

smb_size_anova <- Anova(glm_size_smb_gamma, type = 2) %>%
  as.data.frame()

glm_size_smb <- glmmTMB(
  tl ~ season + location,
  data = df_smb_size,
  family = tweedie(link = "log"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000)
  )
)

car::Anova(glm_size_smb, type = 2)

smb_size_anova <- Anova(glm_size_smb, type = 2) %>%
  as.data.frame()

sim_size_smb <- simulateResiduals(glm_size_smb)
plot(sim_size_smb)

smb_size_contrasts <- as.data.frame(
  pairs(emmeans(glm_size_smb, ~ season + location, type = "response"), adjust = "tukey")
)

df_emm_size_smb <- emmeans(glm_size_smb, ~ season + location, type = "response") %>%
  as.data.frame() %>%
  mutate(
    season = factor(as.character(season), levels = c("Spring", "Summer", "Fall")),
    location = factor(tolower(as.character(location)), levels = c("creek", "lake"))
  )

plot_size_smb <- ggplot(
  df_emm_size_smb,
  aes(x = season, y = response, group = location)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.6,
    color = "#323332",
    alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = location),
    size = 3,
    shape = 21,
    color = "#323332",
    stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name   = "Location",
    values = habitat_colors,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = NULL,
    y = "SMB total length (mm)"
  ) +
  coord_cartesian(clip = "off") +
  base_theme

# =============================================================================
# COMBINE FIGURE
# =============================================================================

figure_5_combined <- ggarrange(
  plot_flux,
  plot_cpue,
  ncol = 1,
  heights = c(1, 1, 1.1, 1.1),
  align = "v"
)

print(figure_5_combined)

figure_6_combined <- ggarrange(
  plot_size_shiner,
  plot_size_smb,
  ncol = 1,
  heights = c(1, 1, 1.1, 1.1),
  align = "v"
)

print(figure_6_combined)

# =============================================================================
# SAVE
# =============================================================================

ggsave(
  file.path(file_path, "Figures/figure_5_flux_cpue_shiner_smb_size.jpg"),
  plot = figure_5_combined,
  width = 6,
  height = 5,
  dpi = 900
)

ggsave(
  file.path(file_path, "Figures/figure_6.jpg"),
  plot = figure_6_combined,
  width = 6,
  height = 5,
  dpi = 900
)
