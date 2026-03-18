# Functions used in this repository

# Author(s): Charlotte Ward
# Version: 2026-03-14

#Calculate Standard Error
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

#Generate emmeans pairwise results for isotopes
print_emm_pairs_iso <- function(model, response_name) {
  cat("\nEstimated marginal means:\n")
  print(emmeans(model, ~ month * location, type = "response"))
  
  cat("\nPairwise comparisons for month within location:\n")
  print(pairs(
    emmeans(model, ~ month | location, type = "response"),
    adjust = "tukey"
  ))
  
  cat("\nPairwise comparisons for location within month:\n")
  print(pairs(
    emmeans(model, ~ location | month, type = "response"),
    adjust = "tukey"
  ))
}

#Generate emmeans pairwise results for Telemetry
print_emm_pairs_tel <- function(model) {
  cat("\nEstimated marginal means:\n")
  print(emmeans(model, ~ Season * release_location, type = "response"))
  
  cat("\nPairwise comparisons for Season within release location:\n")
  print(pairs(
    emmeans(model, ~ Season | release_location, type = "response"),
    adjust = "tukey"
  ))
  
  cat("\nPairwise comparisons for release location within Season:\n")
  print(pairs(
    emmeans(model, ~ release_location | Season, type = "response"),
    adjust = "tukey"
  ))
}


# Calculate lake-derived carbon and trophic position (uses d13C_kilj)
calculate_prop_tp <- function(iso_data, chosen_species, chosen_tissue) {
  months  <- c("may", "august", "october")
  results <- list()
  
  cr_baseline <- iso_data %>% filter(organism == "mayfly", location == "creek")
  lk_baseline <- iso_data %>% filter(organism == "mayfly", location == "lake")
  
  cr_mean_dC <- mean(cr_baseline$d13C_kilj, na.rm = TRUE)
  cr_mean_dN <- mean(cr_baseline$d15N,      na.rm = TRUE)
  lk_mean_dC <- mean(lk_baseline$d13C_kilj, na.rm = TRUE)
  lk_mean_dN <- mean(lk_baseline$d15N,      na.rm = TRUE)
  
  for (month_i in months) {
    df <- iso_data %>%
      filter(
        organism %in% chosen_species,
        tissue == chosen_tissue,
        month  == month_i
      )
    
    lake_carbon <- (df$d13C_kilj - cr_mean_dC) / (lk_mean_dC - cr_mean_dC)
    lake_carbon <- pmin(pmax(lake_carbon, 0.001), 0.999)
    
    trophic_position <- 2 + (
      lake_carbon       * ((df$d15N - cr_mean_dN) / 3.4) +
        (1 - lake_carbon) * ((df$d15N - lk_mean_dN) / 3.4)
    )
    
    df$lake_carbon      <- lake_carbon
    df$trophic_position <- trophic_position
    results[[month_i]]  <- df
  }
  
  bind_rows(results)
}
