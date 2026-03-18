# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)

# ============================================================
# Load and prep data
# ============================================================
dets <- read_csv("/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv") %>%
  mutate(
    location = if_else(rec_ID %in% c("R012", "R013", "R014", "R015") & location == "transition",
                       "creek", location),
    step_time  = case_when(step_number == 1 ~ step1_dur, step_number == 2 ~ step2_dur, step_number == 3 ~ step3_dur),
    time_delay = case_when(step_number == 1 ~ 0, step_number == 2 ~ step1_dur, step_number == 3 ~ step1_dur + step2_dur),
    release_date        = as.Date(release_date),
    tag_on_date         = as.POSIXct(release_date + days(time_delay), tz = "UTC"),
    tag_off_date        = as.POSIXct(release_date + days(time_delay) + days(step_time), tz = "UTC"),
    detection_timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC")
  )

# ============================================================
# False detection filter + drop invalid tag windows
# ============================================================
detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = FALSE) %>%
  filter(passed_filter == 1)

valid_id_times <- detections_filtered %>%
  group_by(id_time) %>%
  summarise(
    first_det = min(detection_timestamp, na.rm = TRUE),
    last_det  = max(detection_timestamp, na.rm = TRUE),
    tag_on    = first(tag_on_date),
    tag_off   = first(tag_off_date),
    .groups   = "drop"
  ) %>%
  filter(!is.na(tag_on), !is.na(tag_off), tag_off > tag_on,
         last_det >= tag_on, first_det <= tag_off) %>%
  pull(id_time)

detections_filtered <- detections_filtered %>% filter(id_time %in% valid_id_times)

# ============================================================
# Build enriched time series
# ============================================================
create_time_series <- function(data, id_time) {
  fish_data  <- data %>% filter(id_time == !!id_time)
  start_time <- fish_data$tag_on_date[1]  + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) return(NULL)
  
  detection_data <- fish_data %>%
    mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%
    distinct(detection_minute, .keep_all = TRUE) %>%
    select(detection_minute, location)
  
  tibble(
    id_time          = id_time,
    transmitter_id   = fish_data$transmitter_id[1],
    date_time        = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = fish_data$release_location[1]
  ) %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    fill(location, .direction = "down")
}

time_series_all <- map(unique(detections_filtered$id_time),
                       ~ create_time_series(detections_filtered, .x)) %>%
  compact() %>%
  bind_rows() %>%
  mutate(Month_Year = format(date_time, "%Y-%m"),
         id_time    = ifelse(id_time == "1576142_3", "1576142_2", id_time)) %>%
  filter(!(id_time == "1576142_2" & Month_Year < "2024-04"))

# ============================================================
# Shared colours
# ============================================================
creek_pink      <- "#F06C57"
lake_blue       <- "#4A90B8"
transition_grey <- "#C8C8C8"

color_scale <- scale_color_manual(
  name   = "Location",
  breaks = c("creek", "lake", "transition"),
  values = c(creek = creek_pink, lake = lake_blue, transition = transition_grey),
  labels = c("Creek", "Lake", "Transition")
)

shared_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title       = element_text(size = 14)
  )

# ============================================================
# Plot 1: RAW detections — all tags
# ============================================================
dets_abacus <- detections_filtered %>%
  filter(!is.na(location)) %>%
  mutate(detection_time = detection_timestamp) %>%
  arrange(transmitter_id, detection_time) %>%
  group_by(transmitter_id) %>%
  mutate(plot_state = if_else(location != lead(location) & !is.na(lead(location)), "transition", location)) %>%
  ungroup() %>%
  mutate(transmitter_id = factor(transmitter_id,
                                 levels = arrange(distinct(., transmitter_id, detection_time) %>%
                                                    group_by(transmitter_id) %>% slice_min(detection_time), detection_time)$transmitter_id))

p_raw <- ggplot(dets_abacus, aes(x = detection_time, y = transmitter_id, color = plot_state)) +
  geom_point(size = 0.7, alpha = 0.8) +
  color_scale +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y",
                   expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Detection date (UTC)", y = "Transmitter ID") +
  shared_theme

dev.new(); print(p_raw)
ggsave("FigureS5_raw.png", plot = p_raw, width = 8, height = 6.5, dpi = 600)

# ============================================================
# Plot 2: ENRICHED time series — all tags
# ============================================================
ts_abacus <- time_series_all %>%
  filter(!is.na(location)) %>%
  arrange(transmitter_id, date_time) %>%
  group_by(transmitter_id) %>%
  mutate(plot_state = if_else(location != lead(location) & !is.na(lead(location)), "transition", location)) %>%
  ungroup() %>%
  mutate(transmitter_id = factor(transmitter_id,
                                 levels = time_series_all %>% group_by(transmitter_id) %>%
                                   summarise(first = min(date_time, na.rm = TRUE), .groups = "drop") %>%
                                   arrange(first) %>% pull(transmitter_id)))

p_enriched <- ggplot(ts_abacus, aes(x = date_time, y = transmitter_id, color = plot_state)) +
  geom_point(size = 0.7, alpha = 0.8) +
  color_scale +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y",
                   expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Date (UTC)", y = "Transmitter ID") +
  shared_theme +
  theme(plot.title = element_text(size = 13, face = "bold"))

dev.new(); print(p_enriched)
ggsave("FigureS6_enriched.png", plot = p_enriched, width = 8, height = 6.5, dpi = 600)
