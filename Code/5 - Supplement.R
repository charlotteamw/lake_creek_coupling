# Supplementary Materials - Abaccus Plots

# Author(s): Charlotte Ward
# Version: 2026-03-16

# ============================================================
# Load Pkgs
# ============================================================
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(dplyr)

file_path <- getwd()

source(file.path(file_path, "/Code/0 - Functions.R"))

# ============================================================
# Load data
# ============================================================
detections_file <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv"
dets <- read_csv(detections_file)

unique(dets$location)

# ============================================================
# Compute step timing, tagging dates, and convert to POSIXct
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
# (catches: step3_dur = 0, wrong release_date year, etc.)
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
    is.na(tag_on)   |
      is.na(tag_off)  |
      tag_off <= tag_on |          # zero-duration window
      last_det  < tag_on  |        # all detections before window starts
      first_det > tag_off          # all detections after window ends
  )

if (nrow(invalid_id_times) > 0) {
  message("\n--- Dropping ", nrow(invalid_id_times), " id_time(s) with invalid/mismatched windows ---")
  invalid_id_times %>%
    select(id_time, tag_on, tag_off, first_det, last_det, step_dur) %>%
    print(n = Inf)
}

valid_id_times <- window_check %>%
  filter(!id_time %in% invalid_id_times$id_time) %>%
  pull(id_time)

detections_filtered <- detections_filtered %>%
  filter(id_time %in% valid_id_times)

# ============================================================
# Define create_time_series() BEFORE calling it
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
    dplyr::mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%
    dplyr::select(detection_minute, location) %>%
    dplyr::distinct(detection_minute, .keep_all = TRUE)
  
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
# === Abacus plot: all tags, RAW detections ===
# ============================================================
creek_blue      <- "#3A74B4"
lake_orange     <- "#D1775E"
transition_grey <- "#C8C8C8"

dets_abacus <- detections_filtered %>%
  filter(!is.na(location)) %>%
  mutate(detection_time = as.POSIXct(detection_timestamp, tz = "UTC")) %>%
  arrange(transmitter_id, detection_time) %>%
  group_by(transmitter_id) %>%
  mutate(
    is_transition = location != lead(location) & !is.na(location) & !is.na(lead(location)),
    plot_state    = case_when(is_transition ~ "transition", TRUE ~ location)
  ) %>%
  ungroup()

tag_levels <- dets_abacus %>%
  group_by(transmitter_id) %>%
  summarise(first_det = suppressWarnings(min(detection_time, na.rm = TRUE)), .groups = "drop") %>%
  arrange(first_det) %>%
  pull(transmitter_id)

dets_abacus <- dets_abacus %>%
  mutate(transmitter_id = factor(transmitter_id, levels = tag_levels))

x_min <- suppressWarnings(min(dets_abacus$detection_time, na.rm = TRUE))
x_max <- suppressWarnings(max(dets_abacus$detection_time, na.rm = TRUE))

p_abacus <- ggplot(dets_abacus, aes(x = detection_time, y = transmitter_id, color = plot_state)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(
    name   = "State",
    breaks = c("creek", "lake", "transition"),
    values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
    labels = c("Creek", "Lake", "Transition")
  ) +
  scale_x_datetime(
    limits      = c(x_min, x_max),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Detection date (UTC)", y = "Tag (transmitter_id)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y     = element_text(size = 14),
    axis.title.x     = element_text(size = 14)
  )

dev.new(); print(p_abacus)

# ============================================================
# === Per-transmitter abacus plots: raw detections by receiver ===
# ============================================================
all_receivers <- paste0("R", sprintf("%03d", 1:30))

dets_receiver <- detections_filtered %>%
  filter(!is.na(rec_ID), !is.na(location)) %>%
  mutate(
    detection_time = as.POSIXct(detection_timestamp, tz = "UTC"),
    receiver       = factor(rec_ID, levels = all_receivers)
  ) %>%
  filter(!is.na(receiver)) %>%
  arrange(transmitter_id, detection_time) %>%
  group_by(transmitter_id) %>%
  mutate(
    is_transition = location != lead(location) & !is.na(location) & !is.na(lead(location)),
    plot_state    = case_when(is_transition ~ "transition", TRUE ~ location)
  ) %>%
  ungroup()

tx_ids   <- unique(as.character(dets_receiver$transmitter_id))
x_min_rx <- min(dets_receiver$detection_time, na.rm = TRUE)
x_max_rx <- max(dets_receiver$detection_time, na.rm = TRUE)

out_dir <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/receiver_abacus_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

walk(tx_ids, function(tid) {
  fish_dets <- dets_receiver %>% filter(as.character(transmitter_id) == tid)
  if (nrow(fish_dets) == 0) return(invisible(NULL))
  
  n_receivers <- n_distinct(fish_dets$receiver)
  
  p <- ggplot(fish_dets, aes(x = detection_time, y = receiver, color = plot_state)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_manual(
      name   = "State",
      breaks = c("creek", "lake", "transition"),
      values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
      labels = c("Creek", "Lake", "Transition")
    ) +
    scale_y_discrete(name = "Receiver", limits = all_receivers, drop = FALSE) +
    scale_x_datetime(
      limits      = c(x_min_rx, x_max_rx),
      date_breaks = "1 month",
      date_labels = "%b %Y",
      expand      = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title    = paste("Transmitter:", tid),
      subtitle = paste("Detected at", n_receivers, "unique receiver(s)"),
      x        = "Detection date (UTC)",
      y        = "Receiver"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y        = element_text(size = 8),
      axis.text.x        = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
      axis.title         = element_text(size = 12),
      plot.title         = element_text(size = 13, face = "bold"),
      plot.subtitle      = element_text(size = 10, color = "grey40")
    )
  
  dev.new(); print(p)
  ggsave(
    file.path(out_dir, paste0("abacus_receiver_", gsub("[^A-Za-z0-9_]", "_", tid), ".png")),
    plot = p, width = 10, height = 6, dpi = 200
  )
})

# ============================================================
# === Abacus plots for ENRICHED time series (time_series_all) ===
# ============================================================
out_dir_enriched <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/enriched_abacus_plots"
if (!dir.exists(out_dir_enriched)) dir.create(out_dir_enriched, recursive = TRUE)

ts_abacus <- time_series_all %>%
  filter(!is.na(location)) %>%
  arrange(id_time, date_time) %>%
  group_by(id_time) %>%
  mutate(
    is_transition = location != lead(location) & !is.na(location) & !is.na(lead(location)),
    plot_state    = case_when(is_transition ~ "transition", TRUE ~ location)
  ) %>%
  ungroup()

id_time_levels <- ts_abacus %>%
  group_by(id_time) %>%
  summarise(first_det = min(date_time, na.rm = TRUE), .groups = "drop") %>%
  arrange(first_det) %>%
  pull(id_time)

ts_abacus <- ts_abacus %>%
  mutate(id_time = factor(id_time, levels = id_time_levels))

x_min_ts <- min(ts_abacus$date_time, na.rm = TRUE)
x_max_ts <- max(ts_abacus$date_time, na.rm = TRUE)

p_abacus_enriched <- ggplot(ts_abacus, aes(x = date_time, y = id_time, color = plot_state)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(
    name   = "State",
    breaks = c("creek", "lake", "transition"),
    values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
    labels = c("Creek", "Lake", "Transition")
  ) +
  scale_x_datetime(
    limits      = c(x_min_ts, x_max_ts),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  labs(title = "All tags — Enriched time series", x = "Date (UTC)", y = "Tag (id_time)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y     = element_text(size = 14),
    axis.title.x     = element_text(size = 14),
    plot.title       = element_text(size = 13, face = "bold")
  )

dev.new(); print(p_abacus_enriched)
ggsave(
  file.path(out_dir_enriched, "abacus_enriched_all_tags.png"),
  plot = p_abacus_enriched, width = 14, height = 8, dpi = 200
)

id_time_ids <- as.character(levels(ts_abacus$id_time))

walk(id_time_ids, function(tid) {
  fish_ts <- ts_abacus %>% filter(as.character(id_time) == tid)
  if (nrow(fish_ts) == 0) return(invisible(NULL))
  
  n_states <- n_distinct(fish_ts$plot_state)
  
  p <- ggplot(fish_ts, aes(x = date_time, y = plot_state, color = plot_state)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_manual(
      name   = "State",
      breaks = c("creek", "lake", "transition"),
      values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
      labels = c("Creek", "Lake", "Transition")
    ) +
    scale_y_discrete(name = "Location", limits = c("creek", "transition", "lake"), drop = FALSE) +
    scale_x_datetime(
      limits      = c(x_min_ts, x_max_ts),
      date_breaks = "1 month",
      date_labels = "%b %Y",
      expand      = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title    = paste("Tag:", tid, "— Enriched time series"),
      subtitle = paste("Observed in", n_states, "unique state(s)"),
      x        = "Date (UTC)",
      y        = "Location"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y        = element_text(size = 10),
      axis.text.x        = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
      axis.title         = element_text(size = 12),
      plot.title         = element_text(size = 13, face = "bold"),
      plot.subtitle      = element_text(size = 10, color = "grey40")
    )
  
  dev.new(); print(p)
  ggsave(
    file.path(out_dir_enriched, paste0("abacus_enriched_", gsub("[^A-Za-z0-9_]", "_", tid), ".png")),
    plot = p, width = 10, height = 4, dpi = 200
  )
})
