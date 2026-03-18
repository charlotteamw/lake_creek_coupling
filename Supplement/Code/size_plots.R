library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(rstatix)   # shapiro_test, t_test, wilcox_test, adjust_pvalue, add_significance
library(ggpubr)    # stat_pvalue_manual
library(tidyr)
library(purrr)

# =========================
# Read + filter
# =========================
df <- read_csv("/Users/charlotteward/Documents/algonquin_minnow/Final_Figures/final_figures_data/fish_catches_all.csv")

df_filtered <- df %>%
  filter(
    species %in% c(
      "smallmouth bass","golden shiner","creek chub","common shiner",
      "pumpkinseed","yellow perch","brook trout","pearl dace",
      "bluntnose minnow","blacknose shiner","northern redbelly dace"
    ),
    !is.na(tl), tl > 0
  ) %>%
  mutate(
    location = str_to_title(location),                 # "lake"/"creek" -> "Lake"/"Creek"
    location = factor(location, levels = c("Lake","Creek"))
  )

cat("\nCounts per species x habitat (raw):\n")
print(df_filtered %>% count(species, location) %>% arrange(species, location))

# =========================
# Remove outliers (1.5 * IQR rule per species × habitat)
# =========================
df_no_outliers <- df_filtered %>%
  group_by(species, location) %>%
  mutate(
    Q1 = quantile(tl, 0.25, na.rm = TRUE),
    Q3 = quantile(tl, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(tl >= lower_bound & tl <= upper_bound) %>%
  ungroup() %>%
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)

message("Outliers removed (IQR per species x habitat): ", nrow(df_filtered) - nrow(df_no_outliers))

cat("\nCounts per species x habitat (outliers removed):\n")
print(df_no_outliers %>% count(species, location) %>% arrange(species, location))

# =========================
# Summary stats WITH and WITHOUT outliers
# =========================
size_summary_with <- df_filtered %>%
  group_by(species, location) %>%
  summarise(
    n = n(),
    mean_tl = round(mean(tl), 1),
    sd_tl   = round(sd(tl), 1),
    min_tl  = min(tl),
    max_tl  = max(tl),
    .groups = "drop"
  )

size_summary_no_out <- df_no_outliers %>%
  group_by(species, location) %>%
  summarise(
    n = n(),
    mean_tl = round(mean(tl), 1),
    sd_tl   = round(sd(tl), 1),
    min_tl  = min(tl),
    max_tl  = max(tl),
    .groups = "drop"
  )

cat("\nSummary WITH outliers:\n"); print(size_summary_with)
cat("\nSummary WITHOUT outliers:\n"); print(size_summary_no_out)

# (Optional) save
# write_csv(size_summary_with,   "size_summary_by_species_habitat_with_outliers.csv")
# write_csv(size_summary_no_out, "size_summary_by_species_habitat_no_outliers.csv")

# =========================
# Habitat comparisons (Lake vs Creek per species)
# Robust to small n and non-normality
# =========================

# Species in BOTH habitats after outlier removal
species_testable <- df_no_outliers %>%
  group_by(species) %>%
  filter(n_distinct(location) == 2) %>%
  pull(species) %>%
  unique()

species_skipped <- setdiff(unique(df_no_outliers$species), species_testable)
if (length(species_skipped)) {
  message("Skipping species (not present in both habitats after outlier removal): ",
          paste(species_skipped, collapse = ", "))
}

df_tests <- df_filtered %>% filter(species %in% species_testable)

# Helper: safe Shapiro per habitat (returns NA p if n<3 or error)
safe_shapiro_by_loc <- function(dat) {
  dat %>%
    group_by(location) %>%
    summarise(
      n  = dplyr::n(),
      p  = tryCatch({
        if (n < 3) NA_real_ else shapiro.test(tl)$p.value
      }, error = function(e) NA_real_),
      .groups = "drop"
    )
}

# Run tests per species with guards
test_results <- df_tests %>%
  group_by(species) %>%
  group_modify(~{
    ns <- .x %>% count(location, name = "n")
    
    # if any group has <3, skip normality and do Wilcoxon
    small_group <- any(ns$n < 3)
    
    if (!small_group) {
      shap <- safe_shapiro_by_loc(.x)
      both_normal <- all(!is.na(shap$p) & shap$p > 0.05)
      
      if (both_normal) {
        out <- tryCatch(
          t_test(.x, tl ~ location, var.equal = FALSE),
          error = function(e) NULL
        )
        method <- "Welch t-test"
      } else {
        out <- tryCatch(
          wilcox_test(.x, tl ~ location, exact = FALSE),
          error = function(e) NULL
        )
        method <- "Wilcoxon rank-sum"
      }
    } else {
      out <- tryCatch(
        wilcox_test(.x, tl ~ location, exact = FALSE),
        error = function(e) NULL
      )
      method <- "Wilcoxon rank-sum (n<3 in a group)"
    }
    
    if (is.null(out) || nrow(out) == 0) {
      tibble(
        .y. = "tl", group1 = "Lake", group2 = "Creek",
        n1 = ns$n[ns$location == "Lake"] %||% NA_integer_,
        n2 = ns$n[ns$location == "Creek"] %||% NA_integer_,
        statistic = NA_real_, p = NA_real_, method = method
      )
    } else {
      out %>%
        mutate(
          method = method,
          n1 = ns$n[ns$location == group1] %||% NA_integer_,
          n2 = ns$n[ns$location == group2] %||% NA_integer_
        )
    }
  }) %>%
  ungroup() %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

cat("\nLake vs Creek tests per species (outliers removed):\n")
print(test_results %>% select(species, method, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif))

# (Optional) save
# write_csv(test_results, "lake_creek_tests_per_species_no_outliers.csv")

# =========================
# Build p-value labels safely
# =========================
ypos <- df_filtered %>%
  group_by(species) %>%
  summarise(y.position = max(tl, na.rm = TRUE) * 1.05, .groups = "drop")

pval_labels <- test_results %>%
  select(species, group1, group2, p, p.adj, p.adj.signif) %>%
  left_join(ypos, by = "species") %>%
  filter(!is.na(p.adj)) %>%
  mutate(y.position = y.position)

# =========================
# Plot with significance annotations
# =========================
lake_creek_colors <- c("Creek" = "#F06C57", "Lake" = "#4A90B8")

ggplot(df_filtered, aes(x = location, y = tl, fill = location)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = lake_creek_colors) +
  theme_bw() +
  labs(
    y = "Total length (mm)",
    fill = "Location",
    x = NULL
  ) +
  ggpubr::stat_pvalue_manual(
    pval_labels,
    label = "p.adj.signif",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    hide.ns = TRUE,
    inherit.aes = FALSE
  )

