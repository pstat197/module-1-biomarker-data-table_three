library(tidyverse)
library(knitr)

# -------------------
# Read headers
# -------------------
var_names <- read_csv('data/biomarker-raw.csv',
                      col_names = FALSE,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, abbreviation = V2) %>%
  na.omit()

# -------------------
# Load data (no trimming); log10 + z-score
# -------------------
biomarker_no_trim <- read_csv('data/biomarker-raw.csv',
                              skip = 2,
                              col_select = -2L,
                              col_names = c('group', 'empty', pull(var_names, abbreviation), 'ados'),
                              na = c('-', '')
) %>%
  filter(!is.na(group)) %>%
  mutate(across(-c(group, ados),
                ~ scale(log10(.x))[, 1])) %>%
  select(group, ados, everything())

# -------------------
# Per-subject outlier count (|z| > 3)
# -------------------
outlier_counts <- biomarker_no_trim %>%
  mutate(outlier_n = rowSums(across(-c(group, ados), ~ abs(.x) > 3), na.rm = TRUE)) %>%
  select(group, ados, outlier_n)

# -------------------
# Group summary (ASD vs TD)
# -------------------
summary_tbl <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    N      = n(),
    Mean   = mean(outlier_n),
    Median = median(outlier_n),
    SD     = sd(outlier_n),
    Min    = min(outlier_n),
    Max    = max(outlier_n),
    .groups = "drop"
  ) %>%
  arrange(group) %>%
  mutate(across(c(Mean, Median, SD), ~ round(.x, 2)))

# -------------------
# Proportions above thresholds
# -------------------
thresholds <- c(10, 25, 50, 100, 200)
props_tbl <- map_dfr(
  thresholds,
  ~ outlier_counts %>%
    group_by(group) %>%
    summarise(
      threshold       = .x,
      num_at_or_above = sum(outlier_n >= .x),
      N               = n(),
      proportion      = num_at_or_above / N,
      .groups = "drop"
    )
) %>%
  arrange(threshold, group) %>%
  mutate(proportion = round(proportion, 3))

# Show tables
print(kable(summary_tbl,
            caption = "Per-subject outlier counts (|z| > 3) after log10 transform and z-scoring (no trimming)"))

print(kable(props_tbl,
            caption = "Proportion of subjects with outlier counts at or above selected thresholds (|z| > 3)"))

# -------------------
# Flag subject-level outliers (Tukey rule per group)
# -------------------
fences <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    Q1  = quantile(outlier_n, 0.25, na.rm = TRUE),
    Q3  = quantile(outlier_n, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    .groups = "drop"
  )

outlier_flagged <- outlier_counts %>%
  left_join(fences, by = "group") %>%
  mutate(is_outlier = outlier_n < lower | outlier_n > upper)

# -------------------
# Boxplot with outliers highlighted
# -------------------
p_box <- ggplot(outlier_flagged, aes(x = group, y = outlier_n, color = is_outlier)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(
    title = "Per-Subject Outlier Counts by Group",
    x = "Group",
    y = "Number of Outliers (|z| > 3)",
    color = "Flagged as outlier participants"
  ) +
  theme_minimal()

print(p_box)

# -------------------
# Histogram by group (faceted)
# -------------------
p_hist <- ggplot(outlier_counts, aes(x = outlier_n, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(
    title = "Histogram of Per-Subject Outlier Counts",
    x = "Number of Outliers (|z| > 3)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_hist)

# -------------------
# Statistical test: Wilcoxon rank-sum
# -------------------
wilcox.test(outlier_n ~ group, data = outlier_counts)


# How many flagged subjects per group?
outlier_totals <- outlier_flagged %>%
  count(group, is_outlier) %>%
  tidyr::complete(group, is_outlier, fill = list(n = 0)) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

print(outlier_totals)

# 2×2 test: is the proportion of flagged subjects different by group?
tab <- table(outlier_flagged$group, outlier_flagged$is_outlier)
print(tab)
chisq.test(tab)        # or fisher.test(tab) if expected counts are small

