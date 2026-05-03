############################################################
# Pooled Shiner Isotope Analysis (Golden + Common)
# Focus: Carbon signature (lake-derived proportion) and TP
############################################################

# Load packages
library(tidyverse)
library(lubridate)
library(FSA)
library(ggpubr)
library(cowplot)

# ---- Load data ----
isotope_data <- read.csv("/Users/charlotteward/Documents/algonquin_minnow/MixSIAR/metadata/iso_metadata.csv", header = TRUE)

# ---- Helper functions ----
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

calculate_prop_tp <- function(iso_data, chosen_species, chosen_tissue) {
  months <- c("may", "august", "october")
  results <- list()
  
  # Baselines
  cr_baseline <- iso_data %>%
    filter(organism %in% c("mayfly"), location == "creek")
  lk_baseline <- iso_data %>%
    filter(organism %in% c("mayfly"), location == "lake")
  
  cr_mean_dC <- mean(cr_baseline$d13C, na.rm = TRUE)
  cr_mean_dN <- mean(cr_baseline$d15N, na.rm = TRUE)
  lk_mean_dC <- mean(lk_baseline$d13C, na.rm = TRUE)
  lk_mean_dN <- mean(lk_baseline$d15N, na.rm = TRUE)
  
  for (month in months) {
    df <- iso_data %>%
      filter(organism %in% chosen_species, tissue == chosen_tissue, month == month)
    
    lake_carbon <- (df$d13C - cr_mean_dC) / (lk_mean_dC - cr_mean_dC)
    lake_carbon <- pmin(pmax(lake_carbon, 0.001), 0.999)
    
    hab_coupling <- 0.5 - abs(0.5 - lake_carbon)
    TP <- 1 + (((lake_carbon * ((df$d15N - cr_mean_dN) / 3.4)) +
                  ((1 - lake_carbon) * ((df$d15N - lk_mean_dN) / 3.4))))
    
    df$lake_carbon <- lake_carbon
    df$hab_coupling <- hab_coupling
    df$trophic_position <- TP
    results[[month]] <- df
  }
  
  bind_rows(results)
}

# ---- Compute metrics for pooled shiners ----
shiner_species <- c("golden shiner", "common shiner")

shiners_liver <- calculate_prop_tp(isotope_data, shiner_species, "liver") %>%
  distinct() %>%
  mutate(month = factor(month, 
                        levels = c("may", "august", "october"),
                        labels = c("Spring", "Summer", "Fall")))

# ---- Summaries (keep both SE and SD; plots will use SD) ----
shiners_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    avg_lake_carbon = mean(lake_carbon, na.rm = TRUE),
    se_lake_carbon  = standard_error(lake_carbon),
    sd_lake_carbon  = sd(lake_carbon, na.rm = TRUE),
    avg_tp          = mean(trophic_position, na.rm = TRUE),
    se_tp           = standard_error(trophic_position),
    sd_tp           = sd(trophic_position, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Plot styling ----
custom_theme <- theme_minimal() +
  theme(
    axis.text.y   = element_text(size = 12),
    axis.text.x   = element_text(size = 12),
    axis.title    = element_text(size = 13),
    plot.title    = element_text(size = 14, face = "bold"),
    axis.line     = element_line(color = "black"),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 12),
    panel.spacing = unit(1.5, "lines"),
    axis.ticks.x  = element_blank()
    # NOTE: legend ON (no legend.position = "none")
  )

# ---- Define colours (fill) ----
fill_scheme <- c("creek" = "#3A74B4",   # yellow
                 "lake"  = "#D1775E")   # blue

# ---- Plot 1: Lake-derived carbon (SE behind points) ----
plot_shiner_carbon <- ggplot(shiners_summary, aes(x = month, y = avg_lake_carbon, group = location)) +
  # line first
  geom_line(aes(), size = 0.6, alpha = 0.4, color = "black") +
  # error bars behind points
  geom_errorbar(aes(ymin = avg_lake_carbon - se_lake_carbon,
                    ymax = avg_lake_carbon + se_lake_carbon),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  # points with outline
  geom_point(aes(fill = location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek","lake"), labels = c("Creek","Lake")) +
  scale_linetype_manual(name = "Location", values = c("creek" = "dashed", "lake" = "dotted"),
                        breaks = c("creek","lake"), labels = c("Creek","Lake")) +
  labs(y = "Proportion Lake-derived Carbon", x = NULL) +
  ylim(0, 1) +
  custom_theme

# ---- Plot 2: Trophic Position (SE behind points) ----
plot_shiner_tp <- ggplot(shiners_summary, aes(x = month, y = avg_tp, group = location)) +
  geom_line(aes(), size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = avg_tp - se_tp,
                    ymax = avg_tp + se_tp),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek","lake"), labels = c("Creek","Lake")) +
  scale_linetype_manual(name = "Location", values = c("creek" = "dashed", "lake" = "dotted"),
                        breaks = c("creek","lake"), labels = c("Creek","Lake")) +
  labs(y = "Trophic Position", x = NULL) +
  ylim(1.0, ) +
  custom_theme

# (Optional) make legend keys look like your points
plot_shiner_carbon <- plot_shiner_carbon +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = "grey20", stroke = 0.6)))

plot_shiner_carbon
plot_shiner_tp <- plot_shiner_tp +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = "grey20", stroke = 0.6)))

# ---- Combine and display with ONE legend ----
combined_plot <- ggarrange(
  plot_shiner_carbon,
  plot_shiner_tp,
  ncol = 1, nrow = 2,
  align = "v",
  heights = c(1, 1),
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot



# ---- Statistics (optional) ----
kruskal.test(lake_carbon ~ month, data = shiners_liver)
kruskal.test(trophic_position ~ month, data = shiners_liver)

dunnTest(lake_carbon ~ month, data = shiners_liver, method = "bonferroni")
dunnTest(trophic_position ~ month, data = shiners_liver, method = "bonferroni")


shiners_liver <- shiners_liver %>% mutate(TL = as.numeric(tl))


lm_carbon_tl <- lm(lake_carbon ~ TL, data = shiners_liver)
summary(lm_carbon_tl)

# prediction grid
tl_seq <- seq(
  min(shiners_liver$TL, na.rm = TRUE),
  max(shiners_liver$TL, na.rm = TRUE),
  length.out = 100
)

pred_df <- data.frame(TL = tl_seq)

# predict on response scale, WITH SE
pred <- predict(lm_carbon_tl, newdata = pred_df, se.fit = TRUE)

pred_df$fit <- pred$fit
pred_df$se  <- pred$se.fit
pred_df$lwr <- pred_df$fit - 1.96 * pred_df$se
pred_df$upr <- pred_df$fit + 1.96 * pred_df$se

# clamp to [0,1] since it's a proportion
pred_df$lwr <- pmax(pred_df$lwr, 0)
pred_df$upr <- pmin(pred_df$upr, 1)

ggplot(shiners_liver, aes(x = TL, y = lake_carbon)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_ribbon(
    data = pred_df,
    aes(x = TL, ymin = lwr, ymax = upr),   # <-- add x = TL here
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_df,
    aes(x = TL, y = fit),
    linewidth = 1.1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Total length (mm)",
    y = "Proportion lake-derived carbon"
  ) +
  theme_minimal()


par(mfrow = c(2,2))
plot(lm_carbon_tl)
par(mfrow = c(1,1))






shapiro.test(residuals(lm_carbon_tl))


anova_mod <- aov(lake_carbon ~ month * location, data = shiners_liver)
summary(anova_mod)


par(mfrow = c(2,2))
plot(anova_mod)
par(mfrow = c(1,1))

shapiro.test(residuals(anova_mod))

library(car)
leveneTest(lake_carbon ~ month * location, data = shiners_liver)

TukeyHSD(anova_mod)

