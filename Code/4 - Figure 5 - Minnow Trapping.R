library(tidyverse)
library(patchwork)
library(grid) 

file_path <- getwd()

source(file.path(file_path, "/Code/0 - Functions.R"))

###### Minnow trap data

mt_data <- read.csv("//minnowtrap_data.csv", header = TRUE) 

# Filter for shiners and clean
shiner_mt_data <- mt_data %>%
  filter(Species %in% c("Golden Shiner", "Common Shiner"), Site != "CR1") %>%
  mutate(
    wt = as.numeric(Wt),
    tl = as.numeric(TL),
    Location = tolower(Location)
  )

# Build all combinations grid
all_sites <- unique(shiner_mt_data$Site)
all_months <- unique(shiner_mt_data$Month)
all_locations <- unique(shiner_mt_data$Location)

site_grid <- expand.grid(
  Site = all_sites,
  Month = all_months,
  Location = all_locations,
  stringsAsFactors = FALSE
)

# Summarize site biomass
site_total <- shiner_mt_data %>%
  group_by(Month, Location, Site) %>%
  summarise(site_biomass = sum(wt, na.rm = TRUE), .groups = "drop")  # use lowercase `wt` here

# Joinand complete
site_total_complete <- site_grid %>%
  left_join(site_total, by = c("Site", "Month", "Location")) %>%
  mutate(
    Month = factor(Month, levels = c("May", "August", "October")),
    season = factor(Month, labels = c("Spring", "Summer", "Fall")),
    Location = factor(Location, levels = c("creek", "lake"))
  ) %>%
  mutate(Location = factor(Location, levels = c("creek", "lake")))

# Summary stats (used in plot)
avg_mt_shiner <- site_total_complete %>%
  group_by(Month, Location) %>%
  summarise(
    avg_biomass = mean(site_biomass, na.rm = TRUE),
    se_biomass = sd(site_biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    season = factor(Month, levels = c("May", "August", "October"),
                    labels = c("Spring", "Summer", "Fall")),
    Location = tolower(Location)
  )



# Custom colors
habitat_colors <- c("creek" = "#F06C57", "lake" = "#4A90B8")

plot_shiner_mt <- ggplot(site_total_complete, aes(
  x = season,
  y = log(site_biomass + 0.01),
  fill = Location,          # keep fill mapped
  group = interaction(season, Location)
  # REMOVED: color = Location
)) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  )+
  geom_boxplot(
    position = position_dodge2(preserve = "single", width = 0.75),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.75,
    size = 0.5,
    color = "#3d3d3d"         # ADD this for black outlines
  ) +
  scale_fill_manual(
    values = habitat_colors,
    name = "Habitat",
    labels = c("Creek", "Lake")
  ) +
  scale_color_manual(
    values = habitat_colors,
    name = "Habitat",
    labels = c("Creek", "Lake")
  ) +
  xlab("Season") +
  ylab("Log Habitat-Specific Biomass (g)") +
  coord_cartesian(clip = "off") +
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
  )+
  coord_cartesian(ylim = c(1.5, 7)) 


plot_shiner_mt




#### directional trap data
trap_data <- read.csv("//directional_traps_all.csv", header = TRUE)

shiner_data <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner"))


shiner_data <- shiner_data %>%
  mutate(
    season_order = factor(season, levels = c("spring", "summer", "fall")),
    wt = as.numeric(wt),
    number = as.numeric(number)  
  )

daily_biomass_shiner <- shiner_data %>%
  group_by(season, day, direction) %>%
  summarize(
    total_daily_biomass = sum(wt, na.rm = TRUE),
    total_daily_abundance = sum(number, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    season_order = factor(season,
                          levels = c("spring", "summer", "fall"),
                          labels = c("Spring", "Summer", "Fall"))
  )


# Season-level summary (mean ± SE for both biomass and abundance)
directional_summary <- daily_biomass_shiner %>%
  group_by(season_order, direction) %>%
  summarise(
    avg_biomass = mean(total_daily_biomass, na.rm = TRUE),
    se_biomass = sd(total_daily_biomass, na.rm = TRUE) / sqrt(n()),
    avg_abundance = mean(total_daily_abundance, na.rm = TRUE),
    se_abundance = sd(total_daily_abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define consistent fill colors
direction_colors <- c("downstream" = "#F06C57", "upstream" = "#eda195")


plot_shiner_biomass <- ggplot(daily_biomass_shiner, aes(
  x = season_order,
  y = total_daily_biomass,
  fill = direction,
  group = interaction(season_order, direction)
)) +
  geom_boxplot(
    aes(color = "#3d3d3d"),
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.45,
    size = 0.5,
    outlier.shape = NA
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = direction_colors,
    name = "Direction",
    labels = c("downstream" = "Downstream", "upstream" = "Upstream")
  ) +
  scale_color_identity() +
  scale_y_log10() +
  xlab("Season") +
  ylab("Total Daily Biomass (log scale)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "#323332"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0.3, 0.3, 1.2, 0.3), "cm"),
    axis.title.x = element_text(vjust = -1.5),
    panel.spacing = unit(1.5, "lines")
  )

plot_shiner_biomass


# Prepare plotting data: reorder direction and set biomass values
directional_summary_plot <- directional_summary %>%
  mutate(
    direction = factor(direction, levels = c("upstream", "downstream")),
    plot_biomass = ifelse(direction == "downstream", -avg_biomass, avg_biomass),
    plot_se = se_biomass
  )

# Fill and outline colors
direction_colors <- c("downstream" = "#636262", "upstream" = "#dbd9d9")

direction_outline_colors <- c("downstream" = "#323332", "upstream" = "#323332")

# Plot
plot_shiner_biomass_bar <- ggplot(directional_summary_plot, aes(
  x = season_order,
  y = plot_biomass,
  fill = direction
)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.6),
    width = 0.6,
    aes(color = direction),
    alpha = 0.8
  ) +
  geom_errorbar(aes(
    ymin = plot_biomass - plot_se,
    ymax = plot_biomass + plot_se,
    color = direction 
  ),
  position = position_dodge(width = 0.6),
  width = 0.2,
  alpha = 1.0,  
  linewidth = 0.6  
  )+
  geom_hline(yintercept = 0, linetype = 1, color = "#3d3d3d", linewidth = 0.5) +
  scale_fill_manual(
    values = direction_colors,
    name = "Direction",
    labels = c("upstream" = "Upstream", "downstream" = "Downstream")
  ) +
  scale_color_manual(
    values = direction_outline_colors,
    guide = "none"
  ) +
  ylab("Daily Biomass Flux in Creek (g)") +
  xlab("Season") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "#3d3d3d"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )+
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  )+
  ylim(-3000,4370)

plot_shiner_biomass_bar

library(scales)

# Custom signed log transform
signed_log10 <- function(x) {
  sign(x) * log10(abs(x) + 1e-3)  # +0.001g for numerical stability
}

# Your exact same plot (ONLY y-scale changed)
plot_shiner_biomass_bar_log<- ggplot(directional_summary_plot, aes(
  x = season_order,
  y = plot_biomass,
  fill = direction
)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.6),
    width = 0.6,
    aes(color = direction),
    alpha = 0.8
  ) +
  geom_errorbar(aes(
    ymin = plot_biomass - plot_se,
    ymax = plot_biomass + plot_se,
    color = direction 
  ),
  position = position_dodge(width = 0.6),
  width = 0.2,
  alpha = 1.0,  
  linewidth = 0.6  
  )+
  geom_hline(yintercept = 0, linetype = 1, color = "#3d3d3d", linewidth = 0.5) +
  scale_fill_manual(
    values = direction_colors,
    name = "Direction",
    labels = c("upstream" = "Upstream", "downstream" = "Downstream")
  ) +
  scale_color_manual(
    values = direction_outline_colors,
    guide = "none"
  ) +
  # REPLACED: scale_y_continuous() → signed log scale
  scale_y_continuous(
    trans = trans_new(
      name = "signed_log10",
      transform = signed_log10,
      inverse = function(x) ifelse(x > 0, 10^x - 0.001, -(10^-x - 0.001))
    ),
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000")
  ) +
  ylab("Daily Biomass Flux in Creek (signed log₁₀ scale)") +
  xlab("Season") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "#3d3d3d"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )+
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) # Removed ylim() to let log scale determine range

plot_shiner_biomass_bar_log

# Remove x-axis from top plot
plot_shiner_biomass_bar <- plot_shiner_biomass_bar +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Add panel border to bottom plot and spacing fix
plot_shiner_mt <- plot_shiner_mt +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x = element_text(vjust = -1.5)
  )


combined_plot <-  plot_shiner_biomass_bar / 
  plot_spacer() / 
 plot_shiner_mt +
  plot_layout(nrow = 3, heights = c(0.9, 0.0000001, 0.9), guides = "collect") +
  plot_annotation(
    tag_levels = "A", 
    theme = theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      legend.title.align = 0,
      plot.margin = margin(t = 12, b = 12)
    )
  )

combined_plot


ggsave(
  filename = "figure_5.png",
  plot = combined_plot,
  width = 7,
  height = 7.5,
  dpi = 600
)


avg_mt_shiner <- site_total_complete %>%
  group_by(Month, Location) %>%
  summarise(
    avg_biomass = mean(site_biomass, na.rm = TRUE),
    se_biomass = sd(site_biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    season = factor(Month, levels = c("May", "August", "October"),
                    labels = c("Spring", "Summer", "Fall")),
    Location = tolower(Location)
  )


library(glue)

avg_mt_shiner %>%
  mutate(summary = glue("{round(avg_biomass, 1)} ± {round(se_biomass, 1)} g")) %>%
  select(season, Location, summary) %>%
  arrange(season, Location)


directional_summary %>%
  mutate(
    direction = factor(direction, levels = c("upstream", "downstream")),
    season = factor(season_order, levels = c("Spring", "Summer", "Fall")),
    biomass_summary = glue("{round(avg_biomass, 1)} ± {round(se_biomass, 1)} g")
  ) %>%
  select(season, direction, biomass_summary) %>%
  arrange(season, direction)



library(car)

# Build two-way ANOVA model
anova_model_A <- aov(log(site_biomass + 1) ~ season * Location, data = site_total_complete)

# 1. Shapiro-Wilk test on residuals
shapiro.test(residuals(anova_model_A))

# 2. QQ plot for residuals
qqnorm(residuals(anova_model_A)); qqline(residuals(anova_model_A))

# 3. Plot residuals vs fitted
plot(anova_model_A, which = 1)

# 4. Levene’s Test for homogeneity of variance
leveneTest(site_biomass ~ interaction(season, Location), data = site_total_complete)



# Build two-way ANOVA model
anova_model_B <- aov(log(total_daily_biomass+1) ~ season * direction, data = daily_biomass_shiner)

anova_model_B

# 1. Shapiro-Wilk test
shapiro.test(residuals(anova_model_B))

# 2. QQ plot
qqnorm(residuals(anova_model_B)); qqline(residuals(anova_model_B))

# 3. Residuals vs fitted
plot(anova_model_B, which = 1)

# 4. Levene’s Test
leveneTest(total_daily_biomass ~ interaction(season, direction), data = daily_biomass_shiner)



summary(anova_model_A)

library(emmeans)

emmeans(anova_model_A, pairwise ~ Location | season)
emmeans(anova_model_A, pairwise ~ season | Location)

# Back-transform emmeans
regrid(emmeans(anova_model_A, ~ Location | season))


summary(anova_model_B)

library(emmeans)

emmeans(anova_model_B, pairwise ~ direction | season)
emmeans(anova_model_B, pairwise ~ season | direction)

# Back-transform emmeans
regrid(emmeans(anova_model_B, ~ direction | season))


# Create combined group
daily_biomass_shiner <- daily_biomass_shiner %>%
  mutate(group = interaction(season_order, direction))

# Kruskal-Wallis test across all groups
kruskal.test(total_daily_biomass ~ group, data = daily_biomass_shiner)

library(FSA)  # for Dunn's test
dunnTest(total_daily_biomass ~ group, data = daily_biomass_shiner, method = "holm")


