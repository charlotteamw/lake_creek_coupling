

### Note - in Chapter 4.4.3 - states 6 candidate models. Here, we have 12, as all models were run with and without a correlation structure. The correlation structure must be included here to account for autocorrelation is residuals. 

# Load libraries
library(sf)
library(glatos)
library(tidyverse)
library(igraph)
library(vegan)
library(mgcv)
library(car)
library(circlize)
library(lubridate)
library(glmmTMB)
library(DHARMa)

detections_file <- "//detections_clean_alldata.csv"
dets <- read_csv(detections_file)

# ------------------------------
# Filter transmitters/receivers
# ------------------------------
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

# ===============================
# Movement variables
# ===============================

daily_movements <- detections_filtered %>%
  group_by(transmitter_id, detection_date, release_location, tl) %>%
  summarise(unique_locations = n_distinct(location), .groups = "drop") %>%
  mutate(
    multiple_locations = ifelse(unique_locations > 1, 1, 0),
    date = as.Date(detection_date),
    day_of_year = yday(date)
  )

# Active tag sequence
daily_movements <- daily_movements %>%
  left_join(
    detections_filtered %>%
      filter(!is.na(tag_on_date) & !is.na(tag_off_date)) %>%
      group_by(transmitter_id) %>%
      summarise(tag_on = min(tag_on_date), tag_off = max(tag_off_date), .groups = "drop") %>%
      mutate(dates = map2(tag_on, tag_off, ~ seq(.x, .y, by = "day"))) %>%
      unnest(dates) %>%
      group_by(dates) %>%
      summarise(active_tags = n(), .groups = "drop") %>%
      rename(date = dates),
    by = "date"
  ) %>%
  filter(active_tags > 2) %>%
  mutate(
    transmitter_id = as.factor(transmitter_id),
    release_location = as.factor(release_location)
  )


###### Models ######

gamm_1 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML")

gamm_2<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id,),
  method = "REML")

gamm_3<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl, 
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_4<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_5<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)



model_list <- list(
  model_1 = gamm_1,
  model_2= gamm_2,
  model_3= gamm_3,
  model_4= gamm_4, 
  model_5= gamm_5
)

# AIC comparison
map_dbl(model_list, function(model) {
  if ("lme" %in% names(model) && !is.null(model$lme)) {
    AIC(model$lme)
  } else {
    AIC(model$gam)
  }
}) %>% sort()


names(model_list) <- paste0("gamm_", 1:5)

# Model formulas
formulas <- c(
  "~ s(day_of_year) + s(tl) + s(release_location) + AR1",
  "~ s(day_of_year) + tl + s(release_location) + AR1",
  "~ s(day_of_year) + tl + AR1",
  "~ s(day_of_year) + s(release_location) + AR1",
  "~ s(day_of_year) + AR1"
)

# Extract metrics
metrics_df <- map_dfr(model_list, function(model) {
  tibble(
    AIC = AIC(model$lme),
    r.squ = summary(model$gam)$r.sq * 100,
    total_edf = sum(summary(model$gam)$s.table[, "edf"])
  )
})

# Build final summary table
model_comparison <- tibble(
  model = names(model_list),
  formula = formulas,
  AIC = round(metrics_df$AIC, 2),
  r.squ = round(metrics_df$r.squ, 2),
  total_edf = round(metrics_df$total_edf, 2)
)

# Save
write.csv(model_comparison, "GAMM_model_comparison.csv")


#### BEST MODEL ####
summary(gamm_5$gam)

###### Plot predictions ######
plot_gamm_predictions <- function(gamm_model, daily_movements, output_path = NULL) {
  new_data <- tibble(
    day_of_year = 1:366,
    tl = median(daily_movements$tl, na.rm = TRUE),
    transmitter_id = levels(daily_movements$transmitter_id)[1],
    release_location = levels(daily_movements$release_location)[1],
    active_tags = rep(1, 366)
  )
  pred <- predict(gamm_model$gam, newdata = new_data, type = "link", se.fit = TRUE)
  new_data <- new_data %>%
    mutate(
      fit_link = pred$fit,
      lower_CI = plogis(fit_link - 1.96 * pred$se.fit),
      upper_CI = plogis(fit_link + 1.96 * pred$se.fit),
      predicted_prob = plogis(fit_link),
      doy_rotated = (day_of_year - 79) %% 365 + 1
    )
  true_doy_labels <- (seq(0, 360, by = 30) + 79 - 1) %% 365 + 1
  p <- ggplot(new_data, aes(x = doy_rotated, y = predicted_prob)) +
    geom_line(size = 1.1, color = "black") +
    geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "gray50") +
    geom_vline(xintercept = c(94, 186, 277), linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(
      name = "Julian Day",
      breaks = seq(0, 360, by = 30),
      labels = true_doy_labels
    ) +
    labs(y = "Prob. of Lake-Creek Detection") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.margin = unit(c(0.3, 0.3, 1.2, 0.3), "cm"),
      axis.title.x = element_text(vjust = -1.5)
    ) +
    ylim(0, 0.8)
  
  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 8.5, height = 4.5, dpi = 300, bg = "transparent")
  } else {
    return(p)
  }
}

plot_gamm <- plot_gamm_predictions(gamm_5, daily_movements)

plot_gamm

plot(gamm_5$gam)

##### Model Diagnostics

k.check(gamm_5$gam)
gam.check(gamm_5$gam)


# Normalized residuals from the LME component (accounts for AR1)
resids_normalized <- residuals(gamm_5$lme, type = "normalized")

par(mfrow = c(1, 2))
acf(resids_normalized,  main = "ACF - Normalized Residuals")
pacf(resids_normalized, main = "pACF - Normalized Residuals")

# Compare to raw residuals (should show more autocorrelation)
resids_raw <- residuals(gamm_5$lme, type = "response")
acf(resids_raw,  main = "ACF - Raw Residuals")
pacf(resids_raw, main = "pACF - Raw Residuals")

concurvity(gamm_5$gam, full = TRUE)   # each term vs. rest of model
concurvity(gamm_5$gam, full = FALSE)  # pairwise between terms


daily_movements$resid <- residuals(gamm_5$gam, type = "deviance")

ggplot(daily_movements, aes(x = transmitter_id, y = resid)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(x = "Transmitter ID", y = "Deviance Residuals",
       title = "Residuals by Individual")

# Full prediction grid across all 366 days
new_data <- tibble(
  day_of_year = 1:366,
  tl = median(daily_movements$tl, na.rm = TRUE),
  transmitter_id = levels(daily_movements$transmitter_id)[1],
  release_location = levels(daily_movements$release_location)[1],
  active_tags = rep(1, 366)
)

# Get predictions
pred <- predict(gamm_5$gam, newdata = new_data, type = "link", se.fit = TRUE)

# Convert to probabilities with CIs
full_pred <- new_data %>%
  mutate(
    fit_link = pred$fit,
    se_link = pred$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    prob       = plogis(fit_link),
    prob_lower = plogis(lower_link),
    prob_upper = plogis(upper_link)
  )

# Find peak (max prob)
peak <- full_pred[which.max(full_pred$prob), ]
cat("Peak probability:\n")
cat(sprintf("Day %d: %.3f (95%% CI: %.3f–%.3f)\n", 
            peak$day_of_year, peak$prob, peak$prob_lower, peak$prob_upper))

# Find trough (min prob)
trough <- full_pred[which.min(full_pred$prob), ]
cat("Minimum probability:\n")
cat(sprintf("Day %d: %.3f (95%% CI: %.3f–%.3f)\n", 
            trough$day_of_year, trough$prob, trough$prob_lower, trough$prob_upper))

# Optional: also find spring/fall peaks (days ~80–140, ~270–320)
spring_peak <- full_pred %>%
  filter(day_of_year >= 80, day_of_year <= 140) %>%
  slice_max(prob, n = 1)

fall_peak <- full_pred %>%
  filter(day_of_year >= 270, day_of_year <= 320) %>%
  slice_max(prob, n = 1)

cat("Spring peak (days 80-140):\n")
cat(sprintf("Day %d: %.3f (95%% CI: %.3f–%.3f)\n", 
            spring_peak$day_of_year, spring_peak$prob, spring_peak$prob_lower, spring_peak$prob_upper))

cat("Fall peak (days 270-320):\n")
cat(sprintf("Day %d: %.3f (95%% CI: %.3f–%.3f)\n", 
            fall_peak$day_of_year, fall_peak$prob, fall_peak$prob_upper, fall_peak$prob_lower))

sessionInfo()
