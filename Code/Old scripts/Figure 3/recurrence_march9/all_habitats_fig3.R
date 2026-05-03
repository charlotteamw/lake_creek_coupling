############################################################
# FIGURE 3: Isotopes (left) + Telemetry (right)
# Left column: lake-derived carbon & trophic position
# Right column: lake residency & recurrence
############################################################

# ==============================
# Load packages
# ==============================
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(FSA)
library(ggpubr)
library(cowplot)

# ==============================
# Helper(s)
# ==============================
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Shared theme + colours
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

fill_scheme <- c("creek" = "#3A74B4", "lake" = "#D1775E")

# =============================================================================
# PART 1: ISOTOPE ANALYSIS (Golden + Common shiners, liver)
# =============================================================================

isotope_data <- read.csv(
  "/Users/charlotteward/Documents/algonquin_minnow/MixSIAR/metadata/iso_metadata.csv",
  header = TRUE
)

calculate_prop_tp <- function(iso_data, chosen_species, chosen_tissue) {
  months  <- c("may", "august", "october")
  results <- list()
  
  cr_baseline <- iso_data %>% filter(organism %in% c("mayfly"), location == "creek")
  lk_baseline <- iso_data %>% filter(organism %in% c("mayfly"), location == "lake")
  
  cr_mean_dC <- mean(cr_baseline$d13C, na.rm = TRUE)
  cr_mean_dN <- mean(cr_baseline$d15N, na.rm = TRUE)
  lk_mean_dC <- mean(lk_baseline$d13C, na.rm = TRUE)
  lk_mean_dN <- mean(lk_baseline$d15N, na.rm = TRUE)
  
  for (month in months) {
    df <- iso_data %>%
      filter(organism %in% chosen_species, tissue == chosen_tissue, month == month)
    
    lake_carbon <- (df$d13C - cr_mean_dC) / (lk_mean_dC - cr_mean_dC)
    lake_carbon <- pmin(pmax(lake_carbon, 0.001), 0.999)
    
    TP <- 1 + (((lake_carbon       * ((df$d15N - cr_mean_dN) / 3.4)) +
                  ((1 - lake_carbon) * ((df$d15N - lk_mean_dN) / 3.4))))
    
    df$lake_carbon      <- lake_carbon
    df$trophic_position <- TP
    results[[month]]    <- df
  }
  bind_rows(results)
}

shiner_species <- c("golden shiner", "common shiner")

shiners_liver <- calculate_prop_tp(isotope_data, shiner_species, "liver") %>%
  distinct() %>%
  mutate(month = factor(month,
                        levels = c("may", "august", "october"),
                        labels = c("Spring", "Summer", "Fall")))

shiners_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    avg_lake_carbon = mean(lake_carbon, na.rm = TRUE),
    se_lake_carbon  = standard_error(lake_carbon),
    avg_tp          = mean(trophic_position, na.rm = TRUE),
    se_tp           = standard_error(trophic_position),
    .groups         = "drop"
  )

# Plot A: Proportion lake-derived carbon
plot_shiner_carbon <- ggplot(shiners_summary,
                             aes(x = month, y = avg_lake_carbon, group = location)) +
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = avg_lake_carbon - se_lake_carbon,
                    ymax = avg_lake_carbon + se_lake_carbon),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme,
                    breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Proportion Lake-derived Carbon", x = NULL) +
  ylim(0, 1) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3,
                                                 color = "grey20", stroke = 0.6)))

# Plot B: Trophic position
plot_shiner_tp <- ggplot(shiners_summary,
                         aes(x = month, y = avg_tp, group = location)) +
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = avg_tp - se_tp, ymax = avg_tp + se_tp),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme,
                    breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Trophic Position", x = NULL) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3,
                                                 color = "grey20", stroke = 0.6)))

# =============================================================================
# PART 2: TELEMETRY ANALYSIS (Residency + Recurrence)
# =============================================================================

detections_file <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv"
dets <- read_csv(detections_file)

transmitters_to_remove <- c("34905", "32999", "42348", "47754", "32993")
recs_to_remove         <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!rec_ID         %in% recs_to_remove)

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

if (nrow(invalid_id_times) > 0) {
  message("\n--- Dropping ", nrow(invalid_id_times),
          " id_time(s) with invalid/mismatched windows ---")
  invalid_id_times %>%
    select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% (window_check %>%
                         filter(!id_time %in% invalid_id_times$id_time) %>%
                         pull(id_time)))

create_time_series <- function(data, id_time) {
  fish_data <- data %>% filter(id_time == !!id_time)
  
  if (nrow(fish_data) == 0 ||
      is.na(fish_data$tag_on_date[1]) ||
      is.na(fish_data$tag_off_date[1])) {
    warning(paste("Skipping id_time:", id_time, "due to missing or invalid dates"))
    return(NULL)
  }
  
  start_time <- fish_data$tag_on_date[1] + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) {
    warning(paste("Skipping id_time:", id_time, "due to invalid time range"))
    return(NULL)
  }
  
  time_series <- tibble(
    id_time          = id_time,
    date_time        = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = fish_data$release_location[1]
  )
  
  detection_data <- fish_data %>%
    mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%
    select(detection_minute, location) %>%
    distinct(detection_minute, .keep_all = TRUE)
  
  time_series %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    tidyr::fill(location, .direction = "down")
}

id_times         <- unique(detections_filtered$id_time)
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

# ------------------------------
# Seasonal residency
# ------------------------------
residency_seasonal <- time_series_all %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_minutes   = sum(!is.na(location)),
    lake_minutes    = sum(location == "lake",  na.rm = TRUE),
    creek_minutes   = sum(location == "creek", na.rm = TRUE),
    lake_residency  = pmin(lake_minutes  / total_minutes, 1),
    creek_residency = pmin(creek_minutes / total_minutes, 1),
    .groups         = "drop"
  )

# ------------------------------
# Recurrence (return to first detected habitat)
# ------------------------------
recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  fish_data <- time_series_all %>%
    filter(id_time == fish_id) %>%
    arrange(date_time)
  
  fish_data$return_event        <- FALSE
  fish_data$recurrence_interval <- NA_real_
  
  # first observed habitat for this fish
  first_loc <- fish_data %>%
    filter(!is.na(location)) %>%
    slice(1) %>%
    pull(location)
  
  # skip fish with no observed location
  if (length(first_loc) == 0 || is.na(first_loc)) {
    recurrence_seasonal_list[[fish_id]] <- fish_data
    next
  }
  
  exited_target  <- FALSE
  last_exit_time <- as.POSIXct(NA, origin = "1970-01-01", tz = "UTC")
  
  for (i in seq_len(nrow(fish_data))) {
    this_loc <- fish_data$location[i]
    this_time <- fish_data$date_time[i]
    
    if (is.na(this_loc)) next
    
    # fish has left its initial habitat
    if (this_loc != first_loc && !exited_target) {
      exited_target  <- TRUE
      last_exit_time <- this_time
    }
    
    # fish has returned to its initial habitat after leaving it
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
    total_returns            = sum(return_event, na.rm = TRUE),
    mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE),
    .groups                  = "drop"
  )

# ------------------------------
# Merge residency + recurrence, drop Winter
# ------------------------------
final_results <- recurrence_seasonal %>%
  left_join(residency_seasonal, by = c("id_time", "Season", "release_location")) %>%
  distinct(id_time, Season, .keep_all = TRUE) %>%
  filter(!is.na(Season)) %>%
  filter(Season != "Winter")

# ------------------------------
# Telemetry summaries (by Season x release_location)
# ------------------------------
telemetry_summary_by_loc <- final_results %>%
  group_by(Season, release_location) %>%
  summarise(
    mean_lake_residency = mean(lake_residency, na.rm = TRUE),
    se_lake_residency   = standard_error(lake_residency),
    mean_recurrence     = mean(total_returns, na.rm = TRUE),
    se_recurrence       = standard_error(total_returns),
    .groups             = "drop"
  ) %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Fall")))

# Plot C: Lake residency
plot_residency <- ggplot(telemetry_summary_by_loc,
                         aes(x = Season, y = mean_lake_residency,
                             group = release_location)) +
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = mean_lake_residency - se_lake_residency,
                    ymax = mean_lake_residency + se_lake_residency),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 3, shape = 21,
             color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme,
                    breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Mean Lake Residency", x = NULL) +
  ylim(0, 1) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3,
                                                 color = "grey20", stroke = 0.6)))

# Plot D: Recurrence (returns to home habitat)
plot_recurrence <- ggplot(telemetry_summary_by_loc,
                          aes(x = Season, y = mean_recurrence,
                              group = release_location)) +
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = mean_recurrence - se_recurrence,
                    ymax = mean_recurrence + se_recurrence),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 3, shape = 21,
             color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Location", values = fill_scheme,
                    breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  labs(y = "Mean Habitat Recurrence", x = NULL) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3,
                                                 color = "grey20", stroke = 0.6)))

plot_recurrence
# =============================================================================
# PART 3: FIGURE 3 – 2x2 GRID (ISOTOPES LEFT, TELEMETRY RIGHT)
# =============================================================================

left_iso <- ggarrange(
  plot_shiner_carbon + theme(legend.position = "none"),
  plot_shiner_tp     + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v", heights = c(1, 1)
)

right_tel <- ggarrange(
  plot_residency  + theme(legend.position = "none"),
  plot_recurrence + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v", heights = c(1, 1)
)

legend_shared <- get_legend(
  plot_shiner_carbon +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = 12),
          legend.text     = element_text(size = 12))
)

top_row <- ggarrange(
  left_iso, right_tel,
  ncol = 2, nrow = 1, widths = c(1, 1)
)

figure_3 <- ggarrange(
  top_row,
  legend_shared,
  ncol = 1, heights = c(10, 1)
)

figure_3



# ------------------------------
# Count recurrence events for each fish
# Return = return to first detected habitat
# ------------------------------
recurrence_count_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  
  fish_data <- time_series_all %>%
    filter(id_time == fish_id) %>%
    arrange(date_time)
  
  # first detected habitat
  first_loc <- fish_data %>%
    filter(!is.na(location)) %>%
    slice(1) %>%
    pull(location)
  
  # skip fish with no non-missing locations
  if (length(first_loc) == 0 || is.na(first_loc)) {
    recurrence_count_list[[fish_id]] <- tibble(
      id_time = fish_id,
      initial_location = NA_character_,
      n_recurrence_events = NA_integer_
    )
    next
  }
  
  exited <- FALSE
  n_returns <- 0L
  
  for (i in seq_len(nrow(fish_data))) {
    this_loc <- fish_data$location[i]
    
    if (is.na(this_loc)) next
    
    # fish leaves initial habitat
    if (this_loc != first_loc && !exited) {
      exited <- TRUE
    }
    
    # fish returns to initial habitat after leaving
    if (this_loc == first_loc && exited) {
      n_returns <- n_returns + 1L
      exited <- FALSE
    }
  }
  
  recurrence_count_list[[fish_id]] <- tibble(
    id_time = fish_id,
    initial_location = first_loc,
    n_recurrence_events = n_returns
  )
}

recurrence_counts_by_fish <- bind_rows(recurrence_count_list)

recurrence_counts_by_fish
