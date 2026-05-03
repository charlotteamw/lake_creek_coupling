# Generalized Additive Models - Shiner Telemetry

# Author(s): Charlotte Ward & Reilly O'Connor
# Version: 2026-03-18

### Note - all models were run with and without a correlation structure to test whether accounting for autocorrelation is residuals is required... 

# Load libraries
library(sf)
library(glatos)
library(tidyverse)
library(igraph)
library(vegan)
library(mgcv)
library(gratia)
library(car)
library(circlize)
library(lubridate)
library(glmmTMB)
library(DHARMa)
library(bbmle)

file_path <- getwd()

source(file.path(file_path, "/Code/0 - Functions.R"))

##### Load and Clean Data #####
detections_file <- file.path(file_path, "Data/detections_clean_alldata.csv")
dets <- read_csv(detections_file)

transmitters_to_remove <- c("34905", "32999")
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

##### Calculate Movement Variables #####
daily_movements <- detections_filtered %>%
  group_by(transmitter_id, detection_date, release_location, tl) %>%
  summarise(unique_locations = n_distinct(location), .groups = "drop") %>%
  mutate(
    multiple_locations = ifelse(unique_locations > 1, 1, 0),
    date = as.Date(detection_date),
    day_of_year = yday(date)
  )

#Active tag sequence
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

str(daily_movements)

###### GAM Models ######
#Test incorporation of autocorrelation term
gamm_0 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "ML")

gamm_1 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "ML")

gamm_2 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  random = list(transmitter_id = ~1),
  method = "ML")

gamm_3 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  random = list(transmitter_id = ~1),
  method = "ML")

AICtab(gamm_0$lme, gamm_1$lme, gamm_2$lme, gamm_3$lme)



library(dplyr)
library(purrr)
library(tibble)

extract_gamm_table <- function(..., model_names = NULL, digits = 3) {
  mods <- list(...)
  
  if (is.null(model_names)) {
    model_names <- as.character(match.call(expand.dots = FALSE)$...)
  }
  
  stopifnot(length(mods) == length(model_names))
  
  out <- map2_dfr(mods, model_names, function(m, nm) {
    
    lme_obj <- m$lme
    gam_obj <- m$gam
    
    ll  <- as.numeric(logLik(lme_obj))
    k   <- attr(logLik(lme_obj), "df")
    n   <- nobs(lme_obj)
    aic <- AIC(lme_obj)
    bic <- BIC(lme_obj)
    
    has_ar1 <- !is.null(lme_obj$modelStruct$corStruct)
    has_re  <- !is.null(lme_obj$modelStruct$reStruct)
    
    tibble(
      model = nm,
      formula = Reduce(paste, deparse(formula(gam_obj))),
      df = k,
      logLik = ll,
      AIC = aic,
      BIC = bic,
      n = n,
      AR1 = has_ar1,
      random_effect = has_re
    )
  })
  
  out <- out %>%
    arrange(AIC) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    ) %>%
    mutate(
      across(where(is.numeric), ~ round(.x, digits))
    )
  
  out
}

comp_autocor <- extract_gamm_table(
  gamm_0, gamm_1, gamm_2, gamm_3,
  model_names = c("gamm_0", "gamm_1", "gamm_2", "gamm_3")
)

comp_autocor
write.csv(comp_autocor, "Results/gamm_autocorrelation_model_comparison.csv", row.names = FALSE)


#Move forward with top model with both AR1 term and random effects
gamm_3a <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  random = list(transmitter_id = ~1),
  method = "ML")

gamm_3b <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl + release_location, 
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  random = list(transmitter_id = ~1),
  method = "ML"
)

gamm_3c <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + release_location,
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  random = list(transmitter_id = ~1),
  method = "ML"
)

gamm_3d <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  random = list(transmitter_id = ~1),
  method = "ML"
)

AICtab(gamm_3a$lme, gamm_3b$lme, gamm_3c$lme, gamm_3d$lme, 
       nobs = nrow(daily_movements), weights = TRUE)

comp_candidates <- extract_gamm_table(
  gamm_3a, gamm_3b, gamm_3c, gamm_3d,
  model_names = c("gamm_3a", "gamm_3b", "gamm_3c", "gamm_3d")
)

comp_candidates
write.csv(comp_candidates, "Results/gamm_candidate_model_comparison.csv", row.names = FALSE)
comp_candidates <- extract_gamm_table(
  gamm_3a, gamm_3b, gamm_3c, gamm_3d,
  model_names = c(
    "s(day_of_year) + s(tl) + release_location",
    "s(day_of_year) + tl + release_location",
    "s(day_of_year) + release_location",
    "s(day_of_year)"
  )
)

comp_candidates


summary(gamm_3c$gam)
summary(gamm_3c$lme)


comp_candidates_pub <- comp_candidates %>%
  transmute(
    Model = model,
    `K` = df,
    `logLik` = logLik,
    `AIC` = AIC,
    `ΔAIC` = delta_AIC,
    Weight = AIC_weight
  )

comp_candidates_pub


#While 3c is a slightly better fit, release_location is non-significant and AIC is within 2, 
#Therefore gamm_3d is top model
summary(gamm_3d$gam)
summary(gamm_3d$lme)
# ranef(gamm_3d$lme)
gam.check(gamm_3d$gam)
draw(gamm_3d$gam)

###### Plot Predicted Lake-Creek Detection Probability ######
new_data <- tibble(
  day_of_year = 1:366,
  # tl = median(daily_movements$tl, na.rm = TRUE),
  transmitter_id = levels(daily_movements$transmitter_id)[1],
  # release_location = levels(daily_movements$release_location)[1],
  active_tags = rep(1, 366)
)

pred <- predict(gamm_3d$gam, newdata = new_data, type = "link", se.fit = TRUE)

new_data <- new_data %>%
  mutate(
    fit_link = pred$fit,
    lower_CI = plogis(fit_link - 1.96 * pred$se.fit),
    upper_CI = plogis(fit_link + 1.96 * pred$se.fit),
    predicted_prob = plogis(fit_link),
    doy_rotated = (day_of_year - 79) %% 365 + 1
  )

true_doy_labels <- (seq(0, 360, by = 30) + 79 - 1) %% 365 + 1

gg_gam <- ggplot(new_data, aes(x = doy_rotated, y = predicted_prob)) +
  geom_line(linewidth = 1.1, color = "black") +
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

gg_gam

ggsave(file.path(file_path, "Figures/figure_4.jpg"), plot = gg_gam, width = 8.5, height = 4.5, dpi = 600, bg = 'transparent')

##### Model Diagnostics #####
k.check(gamm_3d$gam)
gam.check(gamm_3d$gam)

#Look at Normalized residuals from the LME component (accounts for AR1)
resids_normalized <- residuals(gamm_3d$lme, type = "normalized")

#uncorrected
acf(residuals(gamm_3d$lme))
#corrected with AR1
acf(resids_normalized,  main = "ACF - Normalized Residuals")


##### Look at various summary stats, seasonal max min probabilities ####
#Full prediction grid across all 366 days
new_data <- tibble(
  day_of_year = 1:366,
  tl = median(daily_movements$tl, na.rm = TRUE),
  transmitter_id = levels(daily_movements$transmitter_id)[1],
  release_location = levels(daily_movements$release_location)[1],
  active_tags = rep(1, 366)
)

# Get predictions
pred <- predict(gamm_3d$gam, newdata = new_data, type = "link", se.fit = TRUE)

#Convert to probabilities with CIs
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

#Find peak (max prob)
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

