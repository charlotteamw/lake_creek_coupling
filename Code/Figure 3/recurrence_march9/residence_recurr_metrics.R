# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(dplyr)
library(FSA)
library(viridis)

# ============================================================
# Load data
# ============================================================
detections_file <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv"
dets <- read_csv(detections_file)

unique(dets$location)

# ============================================================
# Filter out unwanted transmitters and receivers
# ============================================================
transmitters_to_remove <- c("34905", "32999", "42348")
recs_to_remove <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!rec_ID %in% recs_to_remove)

# ============================================================
# Compute step timing, tagging dates, and convert to POSIXct
# (all in one mutate — avoids column order issues)
# ============================================================
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

# ============================================================
# Apply false detection filter
# ============================================================
detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = FALSE) %>%
  filter(passed_filter == 1)

# ============================================================
# Drop fish whose tag window doesn't overlap their detections
# ============================================================
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
  filter(
    is.na(tag_on)        |
      is.na(tag_off)       |
      tag_off <= tag_on    |
      last_det  < tag_on   |
      first_det > tag_off
  )

if (nrow(invalid_id_times) > 0) {
  message("\n--- Dropping ", nrow(invalid_id_times), " id_time(s) with invalid/mismatched windows ---")
  invalid_id_times %>%
    select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% (window_check %>%
                         filter(!id_time %in% invalid_id_times$id_time) %>%
                         pull(id_time)))

# ============================================================
# Define create_time_series()
# ============================================================
create_time_series <- function(data, id_time) {
  fish_data <- data %>% filter(id_time == !!id_time)
  
  if (nrow(fish_data) == 0 || is.na(fish_data$tag_on_date[1]) || is.na(fish_data$tag_off_date[1])) {
    warning(paste("Skipping id_time:", id_time, "due to missing or invalid dates"))
    return(NULL)
  }
  
  start_time <- fish_data$tag_on_date[1] + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) {
    warning(paste("Skipping id_time:", id_time, "due to invalid time range"))
    return(NULL)
  }
  
  release_loc <- fish_data$release_location[1]
  
  time_series <- tibble(
    id_time          = id_time,
    date_time        = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = release_loc
  )
  
  detection_data <- fish_data %>%
    dplyr::mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%  # ← POSIXct UTC
    dplyr::select(detection_minute, location) %>%
    dplyr::distinct(detection_minute, .keep_all = TRUE)  # ← prevent row explosion
  
  enriched_series <- time_series %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    fill(location, .direction = "down")
  
  return(enriched_series)
}

# ============================================================
# Generate enriched time series from detections_filtered
# ============================================================
id_times <- unique(detections_filtered$id_time)
time_series_list <- map(id_times, ~ create_time_series(detections_filtered, .x))
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

# ============================================================
# Season assignment
# ============================================================
time_series_all <- time_series_all %>%
  mutate(
    Season = case_when(
      (month(date_time) == 3 & day(date_time) >= 21) | month(date_time) == 4 | month(date_time) == 5 | (month(date_time) == 6 & day(date_time) <= 20) ~ "Spring",
      (month(date_time) == 6 & day(date_time) >= 21) | month(date_time) == 7 | month(date_time) == 8 | (month(date_time) == 9 & day(date_time) <= 20) ~ "Summer",
      (month(date_time) == 9 & day(date_time) >= 21) | month(date_time) == 10 | month(date_time) == 11 | (month(date_time) == 12 & day(date_time) <= 20) ~ "Fall",
      (month(date_time) == 12 & day(date_time) >= 21) | month(date_time) %in% c(1, 2) | (month(date_time) == 3 & day(date_time) <= 20) ~ "Winter",
      TRUE ~ NA_character_
    )
  )


# ============================================================
# Compute seasonal residency
# ============================================================
residency_seasonal <- time_series_all %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_minutes   = sum(!is.na(location)),
    lake_minutes    = sum(location == "lake",  na.rm = TRUE),
    creek_minutes   = sum(location == "creek", na.rm = TRUE),
    lake_residency  = pmin(lake_minutes  / total_minutes, 1),
    creek_residency = pmin(creek_minutes / total_minutes, 1),
    .groups = "drop"
  )

# ============================================================
# Recurrence: per-fish return events
# ============================================================
recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  fish_data <- time_series_all %>% filter(id_time == fish_id)
  
  exited <- FALSE
  fish_data$return_event        <- FALSE
  fish_data$recurrence_interval <- NA_real_
  last_exit_time                <- NA
  
  for (i in seq_len(nrow(fish_data))) {
    if (!is.na(fish_data$location[i]) && fish_data$location[i] != "lake" && !exited) {
      exited         <- TRUE
      last_exit_time <- fish_data$date_time[i]
    }
    
    if (!is.na(fish_data$location[i]) && fish_data$location[i] == "lake" && exited) {
      fish_data$return_event[i] <- TRUE
      
      if (!is.na(last_exit_time)) {
        interval <- as.numeric(difftime(fish_data$date_time[i], last_exit_time, units = "mins"))
        fish_data$recurrence_interval[i] <- ifelse(interval > 0, interval, NA)
      }
      
      exited         <- FALSE
      last_exit_time <- NA
    }
  }
  
  recurrence_seasonal_list[[fish_id]] <- fish_data
}

recurrence_seasonal <- bind_rows(recurrence_seasonal_list) %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_returns            = sum(return_event, na.rm = TRUE),
    mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Merge residency + recurrence
# ============================================================
final_results <- recurrence_seasonal %>%
  left_join(residency_seasonal, by = c("id_time", "Season", "release_location")) %>%
  distinct(id_time, Season, .keep_all = TRUE) %>%
  filter(!is.na(Season))

# ============================================================
# Remove Winter from final analyses/plots
# ============================================================
final_results <- final_results %>%
  filter(Season != "Winter")


# ============================================================
# Summary statistics
# ============================================================
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

fish_avg_lake_residency <- final_results %>%
  group_by(Season, release_location) %>%
  reframe(
    avg_lake_residency = mean(lake_residency, na.rm = TRUE),
    se_lake_residency  = standard_error(lake_residency)
  )

fish_avg_creek_residency <- final_results %>%
  group_by(Season, release_location) %>%
  reframe(
    avg_creek_residency = mean(creek_residency, na.rm = TRUE),
    se_creek_residency  = standard_error(creek_residency)
  )

fish_avg_recurr <- final_results %>%
  group_by(Season, release_location) %>%
  reframe(
    avg_recurrence = mean(total_returns, na.rm = TRUE),
    se_recurrence  = standard_error(total_returns)
  )

fish_avg_recurr_int <- final_results %>%
  na.omit() %>%
  group_by(Season, release_location) %>%
  reframe(
    avg_recurrence_interval = mean(mean_recurrence_interval, na.rm = TRUE),
    se_recurrence_interval  = standard_error(mean_recurrence_interval)
  )

fish_avg <- reduce(
  list(fish_avg_lake_residency, fish_avg_creek_residency, fish_avg_recurr, fish_avg_recurr_int),
  left_join,
  by = c("Season", "release_location")
) %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Fall", "Winter")))

telemetry_summary_by_loc <- final_results %>%
  group_by(Season, release_location) %>%
  summarise(
    mean_lake_residency      = mean(lake_residency, na.rm = TRUE),
    se_lake_residency        = standard_error(lake_residency),
    sd_lake_residency        = sd(lake_residency, na.rm = TRUE),
    mean_recurrence          = mean(total_returns, na.rm = TRUE),
    se_recurrence            = standard_error(total_returns),
    sd_recurrence            = sd(total_returns, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Fall", "Winter")))


# ============================================================
# Statistical tests
# ============================================================
kruskal.test(lake_residency  ~ Season, data = final_results)
kruskal.test(creek_residency ~ Season, data = final_results)
kruskal.test(total_returns   ~ Season, data = final_results)

dunnTest(lake_residency  ~ Season, data = final_results, method = "bonferroni")
dunnTest(creek_residency ~ Season, data = final_results, method = "bonferroni")
dunnTest(total_returns   ~ Season, data = final_results, method = "bonferroni")

# ============================================================
# Plots
# ============================================================

scale_color_manual(
  values = c("creek" = "#3A74B4", "lake" = "#D1775E")
)

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

plot_residency <- ggplot(telemetry_summary_by_loc, 
                         aes(x = Season, y = mean_lake_residency, group = release_location)) +
  # line first (thin, faint)
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  # SE error bars behind points
  geom_errorbar(aes(ymin = mean_lake_residency - se_lake_residency,
                    ymax = mean_lake_residency + se_lake_residency),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  # filled points with outline
  geom_point(aes(fill = release_location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Release\nlocation", values = fill_scheme) +
  labs(y = "Mean Lake Residency", x = NULL) +
  ylim(0, 1) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = "grey20", stroke = 0.6)))

plot_recurrence <- ggplot(telemetry_summary_by_loc, 
                          aes(x = Season, y = mean_recurrence, group = release_location)) +
  geom_line(size = 0.6, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = mean_recurrence - se_recurrence,
                    ymax = mean_recurrence + se_recurrence),
                width = 0.15, size = 0.6, color = "black", alpha = 0.7) +
  geom_point(aes(fill = release_location), size = 3, shape = 21, color = "grey20", stroke = 0.6) +
  scale_fill_manual(name = "Release\nlocation", values = fill_scheme) +
  labs(y = "Mean Lake Recurrence", x = NULL) +
  custom_theme +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = "grey20", stroke = 0.6)))


# ---- Combine with ONE legend at bottom ----
combined_telemetry <- ggarrange(
  plot_residency,
  plot_recurrence,
  ncol = 1, nrow = 2,
  align = "v",
  heights = c(1, 1),
  common.legend = TRUE,
  legend = "bottom"
)

combined_telemetry
