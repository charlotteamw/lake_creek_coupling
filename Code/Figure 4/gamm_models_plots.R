

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

# Load data
detections_file <- "/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv"
dets <- read_csv(detections_file)

# Filter out unwanted transmitters and receivers
dets <- dets %>%
  filter(!transmitter_id %in% c("1594949", "1576137"),
         !rec_ID %in% c("R018", "R017", "R016"))

# Compute step timing and tagging dates
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
    tag_on_date = release_date + days(time_delay),
    tag_off_date = tag_on_date + days(step_time),
    detection_date = as.Date(detection_timestamp_utc)
  )

# Apply false detection filter
detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = TRUE) %>%
  filter(passed_filter == 1)

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
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) + s(transmitter_id, bs = "re") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML")

gamm_2 <- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(tl) +  s(transmitter_id, bs = "re") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML")


gamm_3<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl +  s(transmitter_id, bs = "re") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id,),
  method = "REML")

gamm_4<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl +  s(transmitter_id, bs = "re") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML"
)

gamm_5<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl +  s(transmitter_id, bs = "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_6<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + tl +  s(transmitter_id, bs = "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML"
)


gamm_7<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") +  tl,
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_8<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") +  tl,
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML"
)

gamm_9<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_10<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc") + s(release_location, bs= "re"),
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML"
)

gamm_11<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc"),
  family = binomial(link = "logit"),
  data = daily_movements,
  correlation = corAR1(form = ~ day_of_year | transmitter_id),
  method = "REML"
)

gamm_12<- gamm(
  multiple_locations ~ s(day_of_year, bs = "cc"),
  family = binomial(link = "logit"),
  data = daily_movements,
  method = "REML"
)


model_list <- list(
  model_1 = gamm_1,
  model_2= gamm_2,
  model_3= gamm_3,
  model_4= gamm_4, 
  model_5= gamm_5, 
  model_6= gamm_6, 
  model_7= gamm_7, 
  model_8 = gamm_8, 
  model_9 = gamm_9, 
  model_10 = gamm_10,
  model_11 = gamm_11, 
  model_12 = gamm_12
  
)

# AIC comparison
map_dbl(model_list, function(model) {
  if ("lme" %in% names(model) && !is.null(model$lme)) {
    AIC(model$lme)
  } else {
    AIC(model$gam)
  }
}) %>% sort()


names(model_list) <- paste0("gamm_", 1:12)

# Model formulas
formulas <- c(
  "~ s(day_of_year) + s(tl) + s(transmitter_id) + s(release_location) + AR1",
  "~ s(day_of_year) + s(tl) + s(transmitter_id) + s(release_location)",
  "~ s(day_of_year) + tl + s(transmitter_id) + s(release_location) + AR1",
  "~ s(day_of_year) + tl + s(transmitter_id) + s(release_location)",
  "~ s(day_of_year) + tl + s(transmitter_id) + AR1",
  "~ s(day_of_year) + tl + s(transmitter_id)",
  "~ s(day_of_year) + tl + AR1",
  "~ s(day_of_year) + tl",
  "~ s(day_of_year) + s(release_location) + AR1",
  "~ s(day_of_year) + s(release_location)",
  "~ s(day_of_year) + AR1",
  "~ s(day_of_year)"
  
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

# Save to Excel
write.xlsx(model_comparison, "GAMM_model_comparison.xlsx", overwrite = TRUE)


#### BEST MODEL ####
summary(gamm_11$gam)

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

plot_gamm <- plot_gamm_predictions(gamm_11, daily_movements)

plot_gamm

plot(gamm_11$gam)

