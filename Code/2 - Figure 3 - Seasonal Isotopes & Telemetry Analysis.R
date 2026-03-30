# Code to create Figure 3 of Manuscript - Comparing Seasonal Differences in Isotopes and Telemetry

# Author(s): Charlotte Ward & Reilly O'Connor
# Version: 2026-03-16


# Load packages
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

source(file.path(file_path, "/Code/0 - Functions.R"))

# Type III tests require sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

##### Set plot themes & Colours #####
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
  )


fill_scheme <- c("creek" = "#F06C57", "lake" = "#4A90B8")


##### Set Up Isotope Data, summarize mean, se... #####
isotope_data <- read.csv(file.path(file_path, "Data/iso_metadata.csv"), header = T)

shiner_species <- c("golden shiner", "common shiner")

shiners_liver <- calculate_prop_tp(isotope_data, shiner_species, "liver") %>%
  distinct() %>%
  mutate(
    location = factor(location, levels = c("creek", "lake")),
    month    = factor(
      month,
      levels = c("may", "august", "october"),
      labels = c("Spring", "Summer", "Fall")
    )
  )

shiners_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    avg_lake_carbon = mean(lake_carbon,      na.rm = TRUE),
    se_lake_carbon  = standard_error(lake_carbon),
    avg_tp          = mean(trophic_position, na.rm = TRUE),
    se_tp           = standard_error(trophic_position),
    .groups         = "drop"
  )


##### Baselines Summary and Models #####
mayfly_stats <- isotope_data %>%
  filter(organism == "mayfly") %>%
  group_by(location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj,   na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N,      na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N,        na.rm = TRUE), 2),
    .groups   = "drop"
  )
cat("\n--- Mayfly isotopes by location ---\n")
print(mayfly_stats)

mussel_stats <- isotope_data %>%
  filter(organism == "mussel") %>%
  group_by(location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj,   na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N,      na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N,        na.rm = TRUE), 2),
    .groups   = "drop"
  )
cat("\n--- Mussel isotopes by location ---\n")
print(mussel_stats)


#set up data for model, including mayfly and mussel...
baseline_data <- isotope_data %>%
  filter(organism %in% c("mayfly", "mussel")) %>%
  mutate(
    location = factor(location, levels = c("creek", "lake")),
    month = factor(
      month,
      levels = c("may", "august", "october"),
      labels = c("Spring", "Summer", "Fall")
    )
  )

##### Mayfly C13 N15 Models #####
baseline_data_mf <- baseline_data %>% filter(organism == "mayfly")

glm_d13C_mf <- glm(
  d13C_kilj ~ month * location,
  data = baseline_data_mf,
  family = gaussian()
)

car::Anova(glm_d13C_mf, type = 3)
#no interaction, location differences only...

print_emm_pairs_iso(glm_d13C_mf, "d13C_kilj")

glm_d15N_mf <- glm(
  d15N ~ month * location,
  data = baseline_data_mf,
  family = gaussian()
)

car::Anova(glm_d15N_mf, type = 3)

print_emm_pairs_iso(glm_d15N, "d15N")

##### Mussel C13 N15 Models #####
baseline_data_muss <- baseline_data %>% filter(organism == "mussel")

glm_d13C_muss <- glm(
  d13C_kilj ~ month * location,
  data = baseline_data_muss,
  family = gaussian()
)

car::Anova(glm_d13C_muss, type = 3)
print_emm_pairs_iso(glm_d13C_muss, "d13C_kilj")

glm_d15N_muss <- glm(
  d15N ~ month * location,
  data = baseline_data_muss,
  family = gaussian()
)

car::Anova(glm_d15N_muss, type = 3)

print_emm_pairs_iso(glm_d15N_muss, "d15N")



##### Golden Shiner & Common Shiner (Shiner) Data set up #####
shiner_iso_stats <- isotope_data %>%
  filter(
    organism %in% c("golden shiner", "common shiner"),
    tissue == "liver"
  ) %>%
  group_by(organism, location) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C_kilj, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C_kilj,   na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N,      na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N,        na.rm = TRUE), 2),
    .groups   = "drop"
  )

head(shiner_iso_stats)

##### Shiner Raw C13 N15 Models #####
shiners_raw <- isotope_data %>%
  filter(
    organism %in% c("golden shiner", "common shiner"),
    tissue == "liver"
  ) %>%
  mutate(
    location = factor(location, levels = c("creek", "lake")),
    month    = factor(
      month,
      levels = c("may", "august", "october"),
      labels = c("Spring", "Summer", "Fall")
    )
  )

glm_d13C_shnr <- glm(
  d13C_kilj ~ month * location,
  data = shiners_raw,
  family = gaussian()
)

car::Anova(glm_d13C_shnr, type = 3)
print_emm_pairs_iso(glm_d13C_shnr, "d13C_kilj")

glm_d15N_shnr <- glm(
  d15N ~ month * location,
  data = shiners_raw,
  family = gaussian()
)

car::Anova(glm_d15N_shnr, type = 3)
print_emm_pairs_iso(glm_d15N_shnr, "d15N")

##### Shiner Prop Lake Carbon & TP Models #####
#use glmTMB for lake_carbon as proportion (0, 1)
glmb_carbon_shnr <- glmmTMB(
  lake_carbon ~ month * location,
  data = shiners_liver,
  family = beta_family(link = "logit")
)

car::Anova(glmb_carbon_shnr, type = 3)
print_emm_pairs_iso(glmb_carbon_shnr, "lake_carbon")

sim_carb <- simulateResiduals(glmb_carbon_shnr)
plot(sim_carb) 
testOutliers(sim_carb, type = "bootstrap", nBoot = 1000)

glm_tp_shnr <- glm(
  trophic_position ~ month * location,
  data = shiners_liver,
  family = gaussian()
)

car::Anova(glm_tp_shnr, type = 3)
print_emm_pairs_iso(glm_tp_shnr, "trophic_position")

sim_tp <- simulateResiduals(glm_tp_shnr)
plot(sim_tp) 

##### Shiner Summary Stats, Delta Calculations #####
mixing_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    n       = n(),
    mean_lc = round(mean(lake_carbon,      na.rm = TRUE), 3),
    se_lc   = round(standard_error(lake_carbon), 3),
    sd_lc   = round(sd(lake_carbon,        na.rm = TRUE), 3),
    mean_tp = round(mean(trophic_position, na.rm = TRUE), 3),
    se_tp   = round(standard_error(trophic_position), 3),
    sd_tp   = round(sd(trophic_position,   na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("\n--- Mean ± SE: lake carbon & trophic position by season & location ---\n")
print(mixing_summary)

cat("\n--- Delta: season-to-season changes within location (lake carbon) ---\n")
for (loc in c("creek", "lake")) {
  d <- mixing_summary %>% filter(location == loc)
  spr <- d %>% filter(month == "Spring") %>% pull(mean_lc)
  sumr <- d %>% filter(month == "Summer") %>% pull(mean_lc)
  fal <- d %>% filter(month == "Fall")   %>% pull(mean_lc)
  cat(sprintf("\n%s:\n", toupper(loc)))
  cat(sprintf("  Spring -> Summer: Δ = %+.3f\n", sumr - spr))
  cat(sprintf("  Summer -> Fall:   Δ = %+.3f\n", fal - sumr))
  cat(sprintf("  Spring -> Fall:   Δ = %+.3f\n", fal - spr))
}

cat("\n--- Delta: lake vs creek per season (lake carbon) ---\n")
for (seas in c("Spring", "Summer", "Fall")) {
  lk <- mixing_summary %>% filter(month == seas, location == "lake")  %>% pull(mean_lc)
  cr <- mixing_summary %>% filter(month == seas, location == "creek") %>% pull(mean_lc)
  cat(sprintf("  %s - lake minus creek: Δ = %+.3f\n", seas, lk - cr))
}

cat("\n--- Delta: season-to-season changes within location (trophic position) ---\n")
for (loc in c("creek", "lake")) {
  d <- mixing_summary %>% filter(location == loc)
  spr <- d %>% filter(month == "Spring") %>% pull(mean_tp)
  sumr <- d %>% filter(month == "Summer") %>% pull(mean_tp)
  fal <- d %>% filter(month == "Fall")   %>% pull(mean_tp)
  cat(sprintf("\n%s:\n", toupper(loc)))
  cat(sprintf("  Spring -> Summer: Δ = %+.3f\n", sumr - spr))
  cat(sprintf("  Summer -> Fall:   Δ = %+.3f\n", fal - sumr))
  cat(sprintf("  Spring -> Fall:   Δ = %+.3f\n", fal - spr))
}

cat("\n--- Delta: lake vs creek per season (trophic position) ---\n")
for (seas in c("Spring", "Summer", "Fall")) {
  lk <- mixing_summary %>% filter(month == seas, location == "lake")  %>% pull(mean_tp)
  cr <- mixing_summary %>% filter(month == seas, location == "creek") %>% pull(mean_tp)
  cat(sprintf("  %s - lake minus creek: Δ = %+.3f\n", seas, lk - cr))
}

##### Shiner Telemetry Data #####
detections_file <- file.path(file_path, "Data/detections_clean_alldata.csv")
dets <- read_csv(detections_file)

#Remove Dead Fish and receivers at Transition Boundary Edges, 
transmitters_to_remove <- c("34905", "32999", "42348")
recs_to_remove         <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!rec_ID %in% recs_to_remove) %>%
  mutate(
    location = case_when(
      tolower(location) == "transition" ~ "creek",
      TRUE ~ location
    )
  )

dets <- dets %>%
  mutate(
    step_time  = case_when(
      step_number == 1 ~ step1_dur,
      step_number == 2 ~ step2_dur,
      step_number == 3 ~ step3_dur
    ),
    time_delay = case_when(
      step_number == 1 ~ 0,
      step_number == 2 ~ step1_dur,
      step_number == 3 ~ step1_dur + step2_dur
    ),
    release_date        = as.Date(release_date),
    tag_on_date         = release_date + days(time_delay),
    tag_off_date        = tag_on_date  + days(step_time),
    detection_date      = as.Date(detection_timestamp_utc),
    detection_timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
    tag_on_date         = as.POSIXct(tag_on_date,  tz = "UTC"),
    tag_off_date        = as.POSIXct(tag_off_date, tz = "UTC")
  )

detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = FALSE) %>%
  filter(passed_filter == 1)

window_check <- detections_filtered %>%
  group_by(id_time) %>%
  summarise(
    first_det = min(detection_timestamp, na.rm = TRUE),
    last_det  = max(detection_timestamp, na.rm = TRUE),
    tag_on    = first(tag_on_date),
    tag_off   = first(tag_off_date),
    step_dur  = first(step_time),
    .groups   = "drop"
  )


invalid_id_times <- window_check %>%
  filter(is.na(tag_on) | is.na(tag_off) | tag_off <= tag_on |
           last_det < tag_on | first_det > tag_off)

invalid_id_times %>%
  dplyr::select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
  print(n = Inf)


if (nrow(invalid_id_times) > 0) {
  message("\n--- Dropping ", nrow(invalid_id_times),
          " id_time(s) with invalid/mismatched windows ---")
  invalid_id_times %>%
    dplyr::select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% (window_check %>%
                         filter(!id_time %in% invalid_id_times$id_time) %>%
                         pull(id_time)))

id_times <- unique(detections_filtered$id_time)
time_series_list <- purrr::map(id_times, ~ create_time_series(detections_filtered, .x))
names(time_series_list) <- id_times
time_series_list <- time_series_list[!sapply(time_series_list, is.null)]

time_series_all <- bind_rows(time_series_list) %>%
  mutate(
    Month_Year = format(date_time, "%Y-%m"),
    id_time    = ifelse(id_time == "1576142_3", "1576142_2", id_time)
  ) %>%
  filter(!(id_time == "1576142_2" & Month_Year < "2024-04")) %>%
  group_by(id_time) %>%
  filter(date_time >= min(date_time, na.rm = TRUE)) %>%
  ungroup()

time_series_all <- time_series_all %>%
  mutate(
    Season = case_when(
      (month(date_time) == 3 & day(date_time) >= 21) | month(date_time) == 4 |
        month(date_time) == 5 | (month(date_time) == 6 & day(date_time) <= 20) ~ "Spring",
      (month(date_time) == 6 & day(date_time) >= 21) | month(date_time) == 7 |
        month(date_time) == 8 | (month(date_time) == 9 & day(date_time) <= 20) ~ "Summer",
      (month(date_time) == 9 & day(date_time) >= 21) | month(date_time) == 10 |
        month(date_time) == 11 | (month(date_time) == 12 & day(date_time) <= 20) ~ "Fall",
      (month(date_time) == 12 & day(date_time) >= 21) | month(date_time) %in% c(1, 2) |
        (month(date_time) == 3 & day(date_time) <= 20) ~ "Winter",
      TRUE ~ NA_character_
    )
  )

str(time_series_all)

##### Calculate Shiner Residency #####
residency_seasonal <- time_series_all %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_minutes   = sum(!is.na(location)),
    lake_minutes    = sum(location == "lake",  na.rm = TRUE),
    creek_minutes   = sum(location == "creek", na.rm = TRUE),
    lake_residency  = ifelse(total_minutes > 0, lake_minutes / total_minutes, NA_real_),
    creek_residency = ifelse(total_minutes > 0, creek_minutes / total_minutes, NA_real_),
    .groups         = "drop"
  )

##### Calculate Shiner Recurrence #####
recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  fish_data <- time_series_all %>%
    filter(id_time == fish_id) %>%
    arrange(date_time)
  
  fish_data$return_event        <- FALSE
  fish_data$recurrence_interval <- NA_real_
  
  first_loc <- fish_data %>%
    filter(!is.na(location)) %>%
    slice(1) %>%
    pull(location)
  
  if (length(first_loc) == 0 || is.na(first_loc)) {
    recurrence_seasonal_list[[fish_id]] <- fish_data
    next
  }
  
  exited_target  <- FALSE
  last_exit_time <- as.POSIXct(NA, origin = "1970-01-01", tz = "UTC")
  
  for (i in seq_len(nrow(fish_data))) {
    this_loc  <- fish_data$location[i]
    this_time <- fish_data$date_time[i]
    
    if (is.na(this_loc)) next
    
    if (this_loc != first_loc && !exited_target) {
      exited_target  <- TRUE
      last_exit_time <- this_time
    }
    
    if (this_loc == first_loc && exited_target) {
      fish_data$return_event[i] <- TRUE
      interval <- as.numeric(difftime(this_time, last_exit_time, units = "mins"))
      fish_data$recurrence_interval[i] <- ifelse(interval > 0, interval, NA_real_)
      exited_target  <- FALSE
      last_exit_time <- as.POSIXct(NA, origin = "1970-01-01", tz = "UTC")
    }
  }
  
  recurrence_seasonal_list[[fish_id]] <- fish_data
}

recurrence_seasonal <- bind_rows(recurrence_seasonal_list) %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_returns            = sum(return_event,         na.rm = TRUE),
    mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE),
    .groups                  = "drop"
  )

##### Merge Telemetry & General Summary #####
final_results <- recurrence_seasonal %>%
  left_join(
    residency_seasonal,
    by = c("id_time", "Season", "release_location")
  ) %>%
  distinct(id_time, Season, .keep_all = TRUE) %>%
  filter(!is.na(Season), Season != "Winter") %>%
  mutate(
    Season           = factor(Season, levels = c("Spring", "Summer", "Fall")),
    release_location = factor(release_location, levels = c("creek", "lake")),
    season_loc       = paste(Season, release_location, sep = "_")
  ) %>%
  filter(!is.na(lake_minutes), !is.na(total_minutes), total_minutes > 0)

#Create adjusted lake residency to allow for beta regression
n <- final_results$total_minutes
final_results$lake_residency_adj <- (final_results$lake_residency * (n - 1) + 0.5) / n

telemetry_summary <- final_results %>%
  group_by(Season, release_location) %>%
  summarise(
    n                   = n(),
    mean_lake_residency = round(mean(lake_residency, na.rm = TRUE), 3),
    se_lake_residency   = round(standard_error(lake_residency), 3),
    mean_recurrence     = round(mean(total_returns,  na.rm = TRUE), 1),
    se_recurrence       = round(standard_error(total_returns), 1),
    .groups             = "drop"
  )
cat("\n--- Telemetry summary by Season × Release Location ---\n")
print(telemetry_summary)

##### Telemetry Models #####
#Lake Residency - Binomial Family
glm_residency <- glmmTMB(
  lake_residency_adj ~ Season * release_location,
  data   = final_results,
  family = beta_family(link = "logit")
)

summary(glm_residency)
car::Anova(glm_residency, type = 3)
# print_emm_pairs_tel(glm_residency)

sim_residency <- simulateResiduals(glm_residency)
plot(sim_residency) 

#Recurrence - Negative binomial family
glm_recurrence <- glmmTMB(
  total_returns ~ Season * release_location,
  data   = final_results,
  family = nbinom2(link = "log")
)

summary(glm_recurrence)
car::Anova(glm_recurrence, type = 3)
print_emm_pairs_tel(glm_recurrence)

sim_recurrence <- simulateResiduals(glm_recurrence)
plot(sim_recurrence) 

##### Figure 3 — Model estimated marginal means from emmeans() for plotting #####
#3A - Lake Derived Carbon
emm_carbon <- emmeans(glmb_carbon_shnr, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

#3B - Trophic position
emm_tp <- emmeans(glm_tp, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = emmean) %>%
  mutate(lo = est - SE, hi = est + SE)

#3C - Lake residency
emm_residency <- emmeans(glm_residency, ~ Season * release_location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

#3D - Recurrence
emm_recurrence <- emmeans(glm_recurrence, ~ Season * release_location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = est - SE, hi = est + SE)

##### Figure 3 #####
#3A - Lake Derived Carbon
plot_shiner_carbon <- ggplot(
  emm_carbon,
  aes(x = month, y = est, group = location)
) +
  geom_line(linewidth = 0.6, alpha = 0.4, color = "#323332") +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7
  ) +
  geom_point(
    aes(fill = location),
    size = 3, shape = 21, color = "#323332", stroke = 0.6
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Proportion Lake-derived Carbon", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3B - Trophic position
plot_shiner_tp <- ggplot(
  emm_tp,
  aes(x = month, y = est, group = location)
) +
  geom_line(linewidth = 0.6, alpha = 0.4, color = "#323332") +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7
  ) +
  geom_point(
    aes(fill = location),
    size = 3, shape = 21, color = "#323332", stroke = 0.6
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Trophic Position", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3C - Lake residency
plot_residency <- ggplot(
  emm_residency,
  aes(x = Season, y = est, group = release_location)
) +
  geom_line(linewidth = 0.6, alpha = 0.4, color = "#323332") +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7
  ) +
  geom_point(
    aes(fill = release_location),
    size = 3, shape = 21, color = "#323332", stroke = 0.6
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Mean Lake Residency", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3D - Recurrence
plot_recurrence <- ggplot(
  emm_recurrence,
  aes(x = Season, y = est, group = release_location)
) +
  geom_line(linewidth = 0.6, alpha = 0.4, color = "#323332") +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7
  ) +
  geom_point(
    aes(fill = release_location),
    size = 3, shape = 21, color = "#323332", stroke = 0.6
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(y = "Mean Habitat Recurrence", x = NULL) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

left_iso <- ggarrange(
  plot_shiner_carbon + theme(legend.position = "none"),
  plot_shiner_tp     + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

right_tel <- ggarrange(
  plot_residency  + theme(legend.position = "none"),
  plot_recurrence + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

legend_shared <- get_legend(
  plot_shiner_carbon +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    )
)

top_row <- ggarrange(left_iso, right_tel, ncol = 2)

figure_3 <- ggarrange(top_row, legend_shared, ncol = 1, heights = c(10, 1))
figure_3

# ggsave(file.path(file_path, "Figures/figure_3.jpg"), plot = figure_3, width = 8, height = 6.5, dpi = 600)
