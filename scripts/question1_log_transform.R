# scripts/01_log_transform.R
# Q1: Why log-transform protein levels?

INPUT_CSV <- "data/biomarker-raw.csv"  
OUTDIR    <- "results/q1_log"
SAMPLE_N  <- 9                         

# ---------- Packages ----------
req <- c("tidyverse", "e1071")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)
library(tidyverse)
library(e1071)

# ---------- Helpers ----------
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

sanitize_filename <- function(x) {
  x %>% gsub("[^A-Za-z0-9._-]+", "_", ., perl = TRUE) %>% substr(1, 120)
}

log_safe <- function(x) {
  x <- as.numeric(x)
  mn <- suppressWarnings(min(x[is.finite(x)], na.rm = TRUE))
  shift <- if (is.finite(mn) && !is.na(mn) && mn > 0) 0 else (-ifelse(is.finite(mn), mn, 0) + 1e-6)
  log10(x + shift + 1e-9)
}

# ---------- Load & Coerce ----------
raw0 <- readr::read_csv(
  INPUT_CSV,
  na = c("", "NA", "NaN", "null", "NULL"),
  guess_max = 1e6,
  show_col_types = FALSE
)

meta_names <- c(
  "Group","group","Diagnosis","diagnosis","Class","class",
  "ID","Id","id","Subject","subject","Sample","sample",
  "Sex","sex","Gender","gender","Age","age","Batch","batch","Plate","plate",
  "Target Full Name","target full name"
)

# only parse numbers from character columns that are NOT metadata
char_cols <- names(raw0)[vapply(raw0, is.character, logical(1))]
cols_to_parse <- setdiff(char_cols, meta_names)

raw <- raw0 %>%
  mutate(across(all_of(cols_to_parse), readr::parse_number))

# Identify protein columns: numeric, enough non-missing values, and not near-binary
is_good_numeric <- function(v) {
  is.numeric(v) &&
    sum(is.finite(v)) >= max(10, floor(0.5 * length(v))) &&
    dplyr::n_distinct(v[is.finite(v)]) > 5
}
protein_cols <- names(raw)[vapply(raw, is_good_numeric, logical(1))]

message("Detected protein columns: ", length(protein_cols))
stopifnot(length(protein_cols) > 0)

# ---------- Sample for visuals ----------
set.seed(197)
proteins_samp <- head(sample(protein_cols, min(SAMPLE_N, length(protein_cols))), SAMPLE_N)

long <- raw %>%
  select(all_of(proteins_samp)) %>%
  pivot_longer(everything(), names_to = "protein", values_to = "value") %>%
  mutate(log_value = log_safe(value))

# ---------- Plots: Histograms ----------
p_raw <- ggplot(long, aes(value)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free") +
  labs(title = "Raw protein levels", x = "Raw intensity", y = "Count")

p_log <- ggplot(long, aes(log_value)) +
  geom_histogram(bins = 30, na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free") +
  labs(title = "Log10-transformed protein levels", x = "log10(value + shift)", y = "Count")

ggsave(file.path(OUTDIR, "hist_raw_sample.png"), p_raw, width = 12, height = 8, dpi = 200)
ggsave(file.path(OUTDIR, "hist_log_sample.png"), p_log, width = 12, height = 8, dpi = 200)

# ---------- Plots: Q–Q ----------
qq_raw <- ggplot(long, aes(sample = value)) +
  stat_qq(na.rm = TRUE) + stat_qq_line(na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free") +
  labs(title = "Q–Q plot (raw)")

qq_log <- ggplot(long, aes(sample = log_value)) +
  stat_qq(na.rm = TRUE) + stat_qq_line(na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free") +
  labs(title = "Q–Q plot (log10)")

ggsave(file.path(OUTDIR, "qq_raw_sample.png"), qq_raw, width = 12, height = 8, dpi = 200)
ggsave(file.path(OUTDIR, "qq_log_sample.png"), qq_log, width = 12, height = 8, dpi = 200)

# ---------- Skewness per sampled protein ----------
skew_tab_sample <- long %>%
  group_by(protein) %>%
  summarise(
    n = sum(is.finite(value)),
    skew_raw = e1071::skewness(value, na.rm = TRUE, type = 2),
    skew_log = e1071::skewness(log_value, na.rm = TRUE, type = 2),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(skew_raw)))

readr::write_csv(skew_tab_sample, file.path(OUTDIR, "skewness_sample.csv"))
print(skew_tab_sample, n = Inf)

# ---------- Skewness across ALL proteins ----------
long_all <- raw %>%
  select(all_of(protein_cols)) %>%
  pivot_longer(everything(), names_to = "protein", values_to = "value") %>%
  mutate(log_value = log_safe(value))

skew_tab_all <- long_all %>%
  group_by(protein) %>%
  summarise(
    n = sum(is.finite(value)),
    skew_raw = e1071::skewness(value, na.rm = TRUE, type = 2),
    skew_log = e1071::skewness(log_value, na.rm = TRUE, type = 2),
    improved = abs(skew_log) < abs(skew_raw),
    .groups = "drop"
  )

readr::write_csv(skew_tab_all, file.path(OUTDIR, "skewness_all.csv"))

cat("\n--- Overall skewness summary ---\n")
cat("Proteins analyzed:", nrow(skew_tab_all), "\n")
cat("Median skew (raw):", median(skew_tab_all$skew_raw, na.rm = TRUE), "\n")
cat("Median skew (log):", median(skew_tab_all$skew_log, na.rm = TRUE), "\n")
cat("Fraction with improved |skew| after log:",
    mean(skew_tab_all$improved, na.rm = TRUE), "\n")

# ---------- Mean–SD relationship (variance stabilization) ----------
mv_tab <- long_all %>%
  group_by(protein) %>%
  summarise(
    mean_raw = mean(value, na.rm = TRUE),
    sd_raw   = sd(value, na.rm = TRUE),
    mean_log = mean(log_value, na.rm = TRUE),
    sd_log   = sd(log_value, na.rm = TRUE),
    .groups = "drop"
  )

g_raw <- ggplot(mv_tab, aes(mean_raw, sd_raw)) +
  geom_point(alpha = 0.6, na.rm = TRUE) +
  labs(title = "Raw scale: SD vs Mean (per protein)", x = "Mean (raw)", y = "SD (raw)")

g_log <- ggplot(mv_tab, aes(mean_log, sd_log)) +
  geom_point(alpha = 0.6, na.rm = TRUE) +
  labs(title = "Log10 scale: SD vs Mean (per protein)", x = "Mean (log10)", y = "SD (log10)")

ggsave(file.path(OUTDIR, "mean_sd_raw.png"), g_raw, width = 8, height = 6, dpi = 200)
ggsave(file.path(OUTDIR, "mean_sd_log.png"), g_log, width = 8, height = 6, dpi = 200)

# For report (t:-
# Using the raw intensities for 1,317 proteins, the distributions are strongly right-skewed (median skew = 2.03). 
# After applying a log10 transform with a minimal shift for non-positive entries, the distributions become substantially 
# more symmetric ( median skew = 0.23), and approximately 70.5% of proteins show reduced absolute skewness. The SD also 
# increases with the mean on the raw scale but this coupling is attenuated on the log scale, indicating variance 
# stabilization. Because many downstream methods assume approximate normality and homoscedasticity, and because fold-
# changes become additive under a log transform, log-transforming these proteomic measurements is appropriate.
