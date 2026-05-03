# =============================================================================
# Figures 3 & 4 — Seasonal Comparisons: Stable Isotopes & Telemetry
# Authors: Charlotte Ward & Reilly O'Connor
# Updated: 2026-04-14
# Description: Models and plots comparing seasonal patterns in shiner isotope
#              composition (habitat coupling, trophic position) and telemetry-based
#              habitat residency and recurrence, plus separate TL panels.
# =============================================================================

# =============================================================================
# 1. SETUP
# =============================================================================

library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(ggpubr)
library(cowplot)
library(car)
library(emmeans)
library(glmmTMB)
library(DHARMa)

file_path <- getwd()
source(file.path(file_path, "Code/0 - Functions.R"))

options(contrasts = c("contr.sum", "contr.poly"))

ensure_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

custom_theme <- theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 14),
    axis.text.x     = element_text(size = 14),
    axis.text.y     = element_text(size = 12)
  )

fill_scheme <- c("creek" = "#F06C57", "lake" = "#4A90B8")
season_levels <- c("Spring", "Summer", "Fall")
month_levels  <- c("may", "august", "october")
month_labels  <- c("Spring", "Summer", "Fall")

legend_override <- guide_legend(
  override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
)

make_bottom_legend <- function(plot_obj) {
  get_legend(
    plot_obj + theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 12)
    )
  )
}

# keep beta-family responses away from exact 0/1
adjust_beta_response <- function(x) {
  n_obs <- sum(!is.na(x))
  if (n_obs <= 1) return(x)
  (x * (n_obs - 1) + 0.5) / n_obs
}

# =============================================================================
# 2. ISOTOPE DATA — LOAD & PROCESS
# =============================================================================

isotope_data <- read.csv(file.path(file_path, "Data/iso_metadata.csv"), header = TRUE)

shiner_species <- c("golden shiner", "common shiner")

shiners_liver <- calculate_prop_tp(isotope_data, shiner_species, "liver") %>%
  distinct() %>%
  mutate(
    tl               = ensure_numeric(tl),
    location         = factor(location, levels = c("creek", "lake")),
    month            = factor(month, levels = month_levels, labels = month_labels),
    habitat_coupling = (0.5 - abs(lake_carbon - 0.5))/2
  ) %>%
  mutate(
    lake_carbon_adj      = adjust_beta_response(lake_carbon)
  )

shiners_raw <- isotope_data %>%
  filter(organism %in% shiner_species, tissue == "liver") %>%
  mutate(
    tl       = ensure_numeric(tl),
    location = factor(location, levels = c("creek", "lake")),
    month    = factor(month, levels = month_levels, labels = month_labels)
  )

# =============================================================================
# 3. BASELINE ISOTOPES — SUMMARY & MODELS
# =============================================================================

mayfly_stats <- isotope_data %>%
  filter(organism == "mayfly") %>%
  group_by(location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N, na.rm = TRUE), 2),
    .groups = "drop"
  )
cat("\n--- Mayfly isotopes by location ---\n")
print(mayfly_stats)

mussel_stats <- isotope_data %>%
  filter(organism == "mussel") %>%
  group_by(location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N, na.rm = TRUE), 2),
    .groups = "drop"
  )
cat("\n--- Mussel isotopes by location ---\n")
print(mussel_stats)

baseline_data <- isotope_data %>%
  filter(organism %in% c("mayfly", "mussel")) %>%
  mutate(
    location = factor(location, levels = c("creek", "lake")),
    month    = factor(month, levels = month_levels, labels = month_labels)
  )

baseline_mf <- filter(baseline_data, organism == "mayfly")
baseline_muss <- filter(baseline_data, organism == "mussel")

glm_d13C_mf <- glm(d13C_kilj ~ month * location, data = baseline_mf, family = gaussian())
d13C_anova_mf <- as.data.frame(Anova(glm_d13C_mf, type = 3))
d13C_contrasts_mf <- as.data.frame(
  pairs(emmeans(glm_d13C_mf, ~ month * location, type = "response"))
)

glm_d15N_mf <- glm(d15N ~ month * location, data = baseline_mf, family = gaussian())
d15N_anova_mf <- as.data.frame(Anova(glm_d15N_mf, type = 3))
d15N_contrasts_mf <- as.data.frame(
  pairs(emmeans(glm_d15N_mf, ~ month * location, type = "response"))
)

glm_d13C_muss <- glm(d13C_kilj ~ month * location, data = baseline_muss, family = gaussian())
glm_d13C_muss <- update(glm_d13C_muss, . ~ . - month:location)
d13C_anova_muss <- as.data.frame(Anova(glm_d13C_muss, type = 2))

glm_d15N_muss <- glm(d15N ~ month * location, data = baseline_muss, family = gaussian())
glm_d15N_muss <- update(glm_d15N_muss, . ~ . - month:location)
d15N_anova_muss <- as.data.frame(Anova(glm_d15N_muss, type = 2))

# =============================================================================
# 4. SHINER ISOTOPE MODELS
# =============================================================================

shiner_iso_stats <- shiners_raw %>%
  group_by(organism, location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Raw isotopes
glm_d13C_shnr <- glm(d13C_kilj ~ month * location, data = shiners_raw, family = gaussian())
d13C_anova <- as.data.frame(Anova(glm_d13C_shnr, type = 3))
d13C_contrasts <- as.data.frame(
  pairs(emmeans(glm_d13C_shnr, ~ month * location, type = "response"))
)

glm_d15N_shnr <- glm(d15N ~ month * location, data = shiners_raw, family = gaussian())
d15N_anova <- as.data.frame(Anova(glm_d15N_shnr, type = 3))
d15N_contrasts <- as.data.frame(
  pairs(emmeans(glm_d15N_shnr, ~ month * location, type = "response"))
)

# Final isotope response models:
# Habitat coupling is highest when lake and creek carbon are equally represented.
glmb_coupling_shnr <- glmmTMB(
  habitat_coupling_adj ~ month * location,
  data = shiners_liver,
  family = beta_family(link = "logit")
)
coupling_anova <- as.data.frame(Anova(glmb_coupling_shnr, type = 3))
coupling_contrasts <- as.data.frame(
  pairs(emmeans(glmb_coupling_shnr, ~ month * location, type = "response"))
)

sim_coupling <- simulateResiduals(glmb_coupling_shnr)
plot(sim_coupling)
testOutliers(sim_coupling, type = "bootstrap", nBoot = 1000)

glm_tp_shnr <- glm(
  trophic_position ~ month * location,
  data = shiners_liver,
  family = gaussian()
)
tp_anova <- as.data.frame(Anova(glm_tp_shnr, type = 3))
tp_contrasts <- as.data.frame(
  pairs(emmeans(glm_tp_shnr, ~ month * location, type = "response"))
)

sim_tp <- simulateResiduals(glm_tp_shnr)
plot(sim_tp)

mixing_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    n             = n(),
    mean_coupling = round(mean(habitat_coupling, na.rm = TRUE), 3),
    se_coupling   = round(standard_error(habitat_coupling), 3),
    sd_coupling   = round(sd(habitat_coupling, na.rm = TRUE), 3),
    mean_tp       = round(mean(trophic_position, na.rm = TRUE), 3),
    se_tp         = round(standard_error(trophic_position), 3),
    sd_tp         = round(sd(trophic_position, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("\n--- Mean ± SE: habitat coupling & trophic position by season × location ---\n")
print(mixing_summary)

# =============================================================================
# 5. TELEMETRY DATA — LOAD & PROCESS
# =============================================================================

dets <- read_csv(
  file.path(file_path, "Data/detections_clean_alldata.csv"),
  show_col_types = FALSE
)

transmitters_to_remove <- c("34905", "32999", "42348")
recs_to_remove <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(
    !transmitter_id %in% transmitters_to_remove,
    !rec_ID %in% recs_to_remove
  ) %>%
  mutate(
    tl = ensure_numeric(tl),
    location = case_when(
      tolower(location) == "transition" ~ "creek",
      TRUE ~ tolower(location)
    )
  )

dets <- dets %>%
  mutate(
    step_time = case_when(
      step_number == 1 ~ step1_dur,
      step_number == 2 ~ step2_dur,
      step_number == 3 ~ step3_dur
    ),
    time_delay = case_when(
      step_number == 1 ~ 0,
      step_number == 2 ~ step1_dur,
      step_number == 3 ~ step1_dur + step2_dur
    ),
    release_date = as.Date(release_date),
    tag_on_date = as.POSIXct(release_date + days(time_delay), tz = "UTC"),
    tag_off_date = as.POSIXct(release_date + days(time_delay) + days(step_time), tz = "UTC"),
    detection_date = as.Date(detection_timestamp_utc),
    detection_timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC")
  )

tl_lookup <- dets %>%
  distinct(id_time, tl) %>%
  filter(!is.na(tl))

detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = FALSE) %>%
  filter(passed_filter == 1)

window_check <- detections_filtered %>%
  group_by(id_time) %>%
  summarise(
    first_det = min(detection_timestamp, na.rm = TRUE),
    last_det = max(detection_timestamp, na.rm = TRUE),
    tag_on = first(tag_on_date),
    tag_off = first(tag_off_date),
    step_dur = first(step_time),
    .groups = "drop"
  )

invalid_id_times <- window_check %>%
  filter(
    is.na(tag_on) | is.na(tag_off) | tag_off <= tag_on |
      last_det < tag_on | first_det > tag_off
  )

if (nrow(invalid_id_times) > 0) {
  message("\n--- Dropping ", nrow(invalid_id_times), " id_time(s) with invalid/mismatched windows ---")
  invalid_id_times %>%
    dplyr::select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% setdiff(window_check$id_time, invalid_id_times$id_time))

id_times <- unique(detections_filtered$id_time)
time_series_list <- purrr::map(id_times, ~ create_time_series(detections_filtered, .x))
names(time_series_list) <- id_times
time_series_list <- time_series_list[!sapply(time_series_list, is.null)]

time_series_all <- bind_rows(time_series_list) %>%
  mutate(
    Month_Year = format(date_time, "%Y-%m"),
    id_time = ifelse(id_time == "1576142_3", "1576142_2", id_time)
  ) %>%
  filter(!(id_time == "1576142_2" & Month_Year < "2024-04")) %>%
  group_by(id_time) %>%
  filter(date_time >= min(date_time, na.rm = TRUE)) %>%
  ungroup()

time_series_all <- time_series_all %>%
  mutate(
    Season = case_when(
      (month(date_time) == 3 & day(date_time) >= 21) |
        month(date_time) %in% 4:5 |
        (month(date_time) == 6 & day(date_time) <= 20) ~ "Spring",
      (month(date_time) == 6 & day(date_time) >= 21) |
        month(date_time) %in% 7:8 |
        (month(date_time) == 9 & day(date_time) <= 20) ~ "Summer",
      (month(date_time) == 9 & day(date_time) >= 21) |
        month(date_time) %in% 10:11 |
        (month(date_time) == 12 & day(date_time) <= 20) ~ "Fall",
      (month(date_time) == 12 & day(date_time) >= 21) |
        month(date_time) %in% 1:2 |
        (month(date_time) == 3 & day(date_time) <= 20) ~ "Winter",
      TRUE ~ NA_character_
    )
  )

# =============================================================================
# 6. TELEMETRY — RESIDENCY & RECURRENCE
# =============================================================================

release_habitat_lookup <- time_series_all %>%
  group_by(id_time) %>%
  summarise(
    target_habitat = first(na.omit(location)),
    .groups = "drop"
  )

residency_seasonal <- time_series_all %>%
  left_join(release_habitat_lookup, by = "id_time") %>%
  group_by(id_time, Season, target_habitat) %>%
  summarise(
    total_minutes = sum(!is.na(location)),
    habitat_minutes = sum(location == target_habitat, na.rm = TRUE),
    habitat_residency = ifelse(total_minutes > 0, habitat_minutes / total_minutes, NA_real_),
    .groups = "drop"
  ) %>%
  rename(release_location = target_habitat)

recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  fish_data <- time_series_all %>%
    filter(id_time == fish_id) %>%
    arrange(date_time)
  
  fish_data$return_event <- FALSE
  fish_data$recurrence_interval <- NA_real_
  
  target_habitat <- fish_data %>%
    filter(!is.na(location)) %>%
    slice(1) %>%
    pull(location)
  
  if (length(target_habitat) == 0 || is.na(target_habitat)) {
    recurrence_seasonal_list[[fish_id]] <- fish_data
    next
  }
  
  exited_target <- FALSE
  last_exit_time <- as.POSIXct(NA, origin = "1970-01-01", tz = "UTC")
  
  for (i in seq_len(nrow(fish_data))) {
    this_loc <- fish_data$location[i]
    this_time <- fish_data$date_time[i]
    
    if (is.na(this_loc)) next
    
    if (this_loc != target_habitat && !exited_target) {
      exited_target <- TRUE
      last_exit_time <- this_time
    }
    
    if (this_loc == target_habitat && exited_target) {
      interval <- as.numeric(difftime(this_time, last_exit_time, units = "mins"))
      fish_data$return_event[i] <- TRUE
      fish_data$recurrence_interval[i] <- ifelse(interval > 0, interval, NA_real_)
      exited_target <- FALSE
      last_exit_time <- as.POSIXct(NA, origin = "1970-01-01", tz = "UTC")
    }
  }
  
  recurrence_seasonal_list[[fish_id]] <- fish_data
}

recurrence_seasonal <- bind_rows(recurrence_seasonal_list) %>%
  left_join(release_habitat_lookup, by = "id_time") %>%
  group_by(id_time, Season, target_habitat) %>%
  summarise(
    total_returns = sum(return_event, na.rm = TRUE),
    mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(release_location = target_habitat)

# =============================================================================
# 7. TELEMETRY — FINAL DATASET
# =============================================================================

final_results <- recurrence_seasonal %>%
  left_join(residency_seasonal, by = c("id_time", "Season", "release_location")) %>%
  left_join(tl_lookup, by = "id_time") %>%
  distinct(id_time, Season, .keep_all = TRUE) %>%
  filter(
    !is.na(Season), Season != "Winter",
    !is.na(habitat_minutes), !is.na(total_minutes), total_minutes > 0,
    !is.na(tl)
  ) %>%
  mutate(
    tl = ensure_numeric(tl),
    Season = factor(Season, levels = season_levels),
    release_location = factor(release_location, levels = c("creek", "lake")),
    season_loc = paste(Season, release_location, sep = "_"),
    habitat_residency_adj = (habitat_residency * (total_minutes - 1) + 0.5) / total_minutes
  )

cat("\n--- TL coverage in telemetry data ---\n")
cat("Rows missing TL:", sum(is.na(final_results$tl)), "\n")
print(summary(final_results$tl))

telemetry_summary <- final_results %>%
  group_by(Season, release_location) %>%
  summarise(
    n = n(),
    mean_habitat_residency = round(mean(habitat_residency, na.rm = TRUE), 3),
    se_habitat_residency = round(standard_error(habitat_residency), 3),
    mean_recurrence = round(mean(total_returns, na.rm = TRUE), 1),
    se_recurrence = round(standard_error(total_returns), 1),
    .groups = "drop"
  )
cat("\n--- Telemetry summary by Season × Release Habitat ---\n")
print(telemetry_summary)

# =============================================================================
# 8. TL MODELS — TL AS RESPONSE
# =============================================================================

glm_tl_iso <- glm(tl ~ month * location, data = shiners_liver, family = gaussian())
tl_iso_anova <- as.data.frame(Anova(glm_tl_iso, type = 3))
tl_iso_contrasts <- as.data.frame(
  pairs(emmeans(glm_tl_iso, ~ month * location, type = "response"), adjust = "tukey")
)

emm_tl_iso <- emmeans(glm_tl_iso, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = emmean) %>%
  mutate(lo = est - SE, hi = est + SE)

glm_tl_tel <- glm(tl ~ Season + release_location, data = final_results, family = gaussian())
tl_tel_anova <- as.data.frame(Anova(glm_tl_tel, type = 2))
tl_tel_contrasts <- as.data.frame(
  pairs(emmeans(glm_tl_tel, ~ Season * release_location, type = "response"), adjust = "tukey")
)

emm_tl_tel <- emmeans(glm_tl_tel, ~ Season * release_location, type = "response") %>%
  as.data.frame() %>%
  rename(est = emmean) %>%
  mutate(lo = est - SE, hi = est + SE)

# =============================================================================
# 9. TELEMETRY RESPONSE MODELS
# =============================================================================

glm_residency <- glmmTMB(
  habitat_residency_adj ~ Season * release_location,
  data = final_results,
  family = beta_family(link = "logit")
)
car::Anova(glm_residency, type = 3)
glm_residency <- update(glm_residency, . ~ . - Season:release_location)
residency_anova <- as.data.frame(Anova(glm_residency, type = 2))

sim_residency <- simulateResiduals(glm_residency)
plot(sim_residency)

glm_recurrence <- glmmTMB(
  total_returns ~ Season * release_location,
  data = final_results,
  family = nbinom2(link = "log")
)
car::Anova(glm_recurrence, type = 3)
glm_recurrence <- update(glm_recurrence, . ~ . - Season:release_location)
recurrence_anova <- as.data.frame(Anova(glm_recurrence, type = 2))

sim_recurrence <- simulateResiduals(glm_recurrence)
plot(sim_recurrence)

# =============================================================================
# 10. ESTIMATED MARGINAL MEANS
# =============================================================================

emm_coupling <- emmeans(glmb_coupling_shnr, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

emm_tp <- emmeans(glm_tp_shnr, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = emmean) %>%
  mutate(lo = est - SE, hi = est + SE)

emm_residency <- emmeans(glm_residency, ~ Season * release_location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

emm_recurrence <- emmeans(glm_recurrence, ~ Season * release_location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = est - SE, hi = est + SE)

# =============================================================================
# 11. TL PANELS FOR FIGURE 4
# =============================================================================
# Keeping Figure 4 as TL-based panels using the original isotope metric (lake carbon)
# unless you also want coupling substituted there.

glmb_carbon_tl <- glmmTMB(
  lake_carbon_adj ~ tl,
  data = shiners_liver,
  family = beta_family(link = "logit")
)

emm_carbon_tl <- tibble(
  tl = seq(min(shiners_liver$tl, na.rm = TRUE), max(shiners_liver$tl, na.rm = TRUE), length.out = 100)
)

pred_carbon_tl <- predict(glmb_carbon_tl, newdata = emm_carbon_tl, type = "link", se.fit = TRUE)
emm_carbon_tl <- emm_carbon_tl %>%
  mutate(
    est = plogis(pred_carbon_tl$fit),
    lo  = plogis(pred_carbon_tl$fit - 1.96 * pred_carbon_tl$se.fit),
    hi  = plogis(pred_carbon_tl$fit + 1.96 * pred_carbon_tl$se.fit)
  )

glm_tp_tl <- glm(trophic_position ~ tl, data = shiners_liver, family = gaussian())

emm_tp_tl <- tibble(
  tl = seq(min(shiners_liver$tl, na.rm = TRUE), max(shiners_liver$tl, na.rm = TRUE), length.out = 100)
)

pred_tp_tl <- predict(glm_tp_tl, newdata = emm_tp_tl, type = "link", se.fit = TRUE)
emm_tp_tl <- emm_tp_tl %>%
  mutate(
    est = pred_tp_tl$fit,
    lo  = pred_tp_tl$fit - 1.96 * pred_tp_tl$se.fit,
    hi  = pred_tp_tl$fit + 1.96 * pred_tp_tl$se.fit
  )

glmb_residency_tl <- glmmTMB(
  habitat_residency_adj ~ tl,
  data = final_results,
  family = beta_family(link = "logit")
)

emm_residency_tl <- tibble(
  tl = seq(min(final_results$tl, na.rm = TRUE), max(final_results$tl, na.rm = TRUE), length.out = 100)
)

pred_residency_tl <- predict(glmb_residency_tl, newdata = emm_residency_tl, type = "link", se.fit = TRUE)
emm_residency_tl <- emm_residency_tl %>%
  mutate(
    est = plogis(pred_residency_tl$fit),
    lo  = plogis(pred_residency_tl$fit - 1.96 * pred_residency_tl$se.fit),
    hi  = plogis(pred_residency_tl$fit + 1.96 * pred_residency_tl$se.fit)
  )

# =============================================================================
# 12. BUILD PANELS
# =============================================================================

# Figure 3: seasonal panels only
plot_shiner_coupling <- ggplot(emm_coupling, aes(x = month, y = est, group = location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Habitat Coupling", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 0.5)) +
  custom_theme + guides(fill = legend_override)

plot_shiner_tp <- ggplot(emm_tp, aes(x = month, y = est, group = location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Trophic Position", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme + guides(fill = legend_override)

plot_recurrence <- ggplot(emm_recurrence, aes(x = Season, y = est, group = release_location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Habitat Recurrence", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme + guides(fill = legend_override)

plot_residency <- ggplot(emm_residency, aes(x = Season, y = est, group = release_location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Habitat Residency", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme + guides(fill = legend_override)

# Figure 4: TL panels only
plot_carbon_tl <- ggplot() +
  geom_ribbon(data = emm_carbon_tl, aes(x = tl, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = emm_carbon_tl, aes(x = tl, y = est), linewidth = 0.8) +
  geom_point(data = shiners_liver, aes(x = tl, y = lake_carbon, fill = location),
             size = 2, shape = 21, color = "#323332", stroke = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Proportion Lake Carbon", x = "Total Length (mm)") +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme

plot_tp_tl <- ggplot() +
  geom_ribbon(data = emm_tp_tl, aes(x = tl, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = emm_tp_tl, aes(x = tl, y = est), linewidth = 0.8) +
  geom_point(data = shiners_liver, aes(x = tl, y = trophic_position, fill = location),
             size = 2, shape = 21, color = "#323332", stroke = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Trophic Position", x = "Total Length (mm)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme

plot_residency_tl <- ggplot() +
  geom_ribbon(data = emm_residency_tl, aes(x = tl, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = emm_residency_tl, aes(x = tl, y = est), linewidth = 0.8) +
  geom_point(data = final_results, aes(x = tl, y = habitat_residency, fill = release_location),
             size = 2, shape = 21, color = "#323332", stroke = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Habitat Residency", x = "Total Length (mm)") +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme

plot_recurrence_tl <- ggplot() +
  geom_point(data = final_results, aes(x = tl, y = total_returns, fill = release_location),
             size = 2, shape = 21, color = "#323332", stroke = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Habitat Recurrence", x = "Total Length (mm)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme

# Optional TL seasonal EMM panels
plot_tl_iso <- ggplot(emm_tl_iso, aes(x = month, y = est, group = location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Total Length (mm)", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme

plot_tl_tel <- ggplot(emm_tl_tel, aes(x = Season, y = est, group = release_location)) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 2.8, shape = 21, color = "#323332", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme, breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Total Length (mm)", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme

# =============================================================================
# 13. ASSEMBLE & SAVE FIGURES
# =============================================================================

fig3_left <- ggarrange(
  plot_shiner_coupling + theme(legend.position = "none"),
  plot_shiner_tp + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

fig3_right <- ggarrange(
  plot_recurrence + theme(legend.position = "none"),
  plot_residency + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

figure_3 <- ggarrange(
  ggarrange(fig3_left, fig3_right, ncol = 2),
  make_bottom_legend(plot_shiner_coupling),
  ncol = 1,
  heights = c(10, 1)
)

figure_3
ggsave(
  file.path(file_path, "Figures/figure_3.jpg"),
  plot = figure_3,
  width = 8,
  height = 7,
  dpi = 900
)

fig4_left <- ggarrange(
  plot_carbon_tl + theme(legend.position = "none"),
  plot_tp_tl + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

fig4_right <- ggarrange(
  plot_residency_tl + theme(legend.position = "none"),
  plot_recurrence_tl + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

figure_4 <- ggarrange(
  ggarrange(fig4_left, fig4_right, ncol = 2),
  make_bottom_legend(plot_carbon_tl),
  ncol = 1,
  heights = c(10, 1)
)

figure_4
ggsave(
  file.path(file_path, "Figures/figure_4.jpg"),
  plot = figure_4,
  width = 8,
  height = 7,
  dpi = 900
)

figure_tl <- ggarrange(
  plot_tl_iso + theme(legend.position = "none"),
  plot_tl_tel + theme(legend.position = "none"),
  make_bottom_legend(plot_tl_iso),
  ncol = 2,
  nrow = 2,
  heights = c(4, 1)
)

figure_tl
ggsave(
  file.path(file_path, "Figures/figure_tl_emms.jpg"),
  plot = figure_tl,
  width = 9,
  height = 5,
  dpi = 900
)