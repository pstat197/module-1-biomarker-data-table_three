# Automatically install missing packages
required_packages <- c(
    "tidyverse", "infer", "randomForest", "tidymodels",
    "modelr", "yardstick", "glmnet"
)

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
} else {
    cat("All required packages are already installed.\n")
}

# Load libraries
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet) # for LASSO

load("data/biomarker-clean.RData")

# Set seed for reproducibility
set.seed(101422)

cat(strrep("=", 70), "\n")
cat("TASK 4: Finding a Simpler Panel with Comparable Accuracy\n")
cat("Method: LASSO Regularization\n")
cat(strrep("=", 70), "\n\n")

# ===============================================================================
# STEP 1: BASELINE ANALYSIS (Replicate in-class approach)
# ===============================================================================

cat("STEP 1: Establishing Baseline Performance\n")
cat(strrep("-", 70), "\n")

## Multiple Testing (T-tests)
test_fn <- function(.df) {
    t_test(.df,
        formula = level ~ group,
        order = c("ASD", "TD"),
        alternative = "two-sided",
        var.equal = F
    )
}

ttests_out <- biomarker_clean %>%
    select(-ados) %>%
    pivot_longer(-group,
        names_to = "protein",
        values_to = "level"
    ) %>%
    nest(data = c(level, group)) %>%
    mutate(ttest = map(data, test_fn)) %>%
    unnest(ttest) %>%
    arrange(p_value) %>%
    mutate(
        m = n(),
        hm = log(m) + 1 / (2 * m) - digamma(1),
        rank = row_number(),
        p.adj = m * hm * p_value / rank
    )

proteins_s1 <- ttests_out %>%
    slice_min(p.adj, n = 10) %>%
    pull(protein)

cat("Top 10 proteins from t-tests:\n")
print(proteins_s1)
cat("\n")

## Random Forest
predictors <- biomarker_clean %>%
    select(-c(group, ados))

response <- biomarker_clean %>%
    pull(group) %>%
    factor()

set.seed(101422)
rf_out <- randomForest(
    x = predictors,
    y = response,
    ntree = 1000,
    importance = T
)

proteins_s2 <- rf_out$importance %>%
    as_tibble() %>%
    mutate(protein = rownames(rf_out$importance)) %>%
    slice_max(MeanDecreaseGini, n = 10) %>%
    pull(protein)

cat("Top 10 proteins from Random Forest:\n")
print(proteins_s2)
cat("\n")

## Baseline: Intersection approach
proteins_baseline <- intersect(proteins_s1, proteins_s2)

cat("Baseline panel (intersection):\n")
cat("Number of proteins:", length(proteins_baseline), "\n")
print(proteins_baseline)
cat("\n")

# Fit baseline model
biomarker_baseline <- biomarker_clean %>%
    select(group, any_of(proteins_baseline)) %>%
    mutate(class = (group == "ASD")) %>%
    select(-group)

set.seed(101422)
biomarker_split_baseline <- biomarker_baseline %>%
    initial_split(prop = 0.8)

fit_baseline <- glm(class ~ .,
    data = training(biomarker_split_baseline),
    family = "binomial"
)

class_metrics <- metric_set(
    sensitivity,
    specificity,
    accuracy,
    roc_auc
)

baseline_results <- testing(biomarker_split_baseline) %>%
    add_predictions(fit_baseline, type = "response") %>%
    mutate(
        pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE)),
        truth_class = factor(class, levels = c(FALSE, TRUE))
    ) %>%
    class_metrics(
        estimate = pred_class,
        truth = truth_class,
        pred,
        event_level = "second"
    )

cat("Baseline Performance Metrics:\n")
print(baseline_results)
cat("\n\n")

# ===============================================================================
# STEP 2: EXPANDED CANDIDATE SET (Union instead of Intersection)
# ===============================================================================

cat("STEP 2: Expanding Candidate Protein Set\n")
cat(strrep("-", 70), "\n")

# Use union to get more candidates for LASSO to select from
proteins_union <- union(proteins_s1, proteins_s2)

cat("Candidate proteins (union of top 10 from each method):\n")
cat("Number of proteins:", length(proteins_union), "\n")
print(proteins_union)
cat("\n\n")

# ===============================================================================
# STEP 3: LASSO REGULARIZATION FOR FEATURE SELECTION
# ===============================================================================

cat("STEP 3: Applying LASSO Regularization\n")
cat(strrep("-", 70), "\n")

# Prepare data with union of proteins
biomarker_union <- biomarker_clean %>%
    select(group, any_of(proteins_union)) %>%
    mutate(class = (group == "ASD")) %>%
    select(-group)

# Use the same train/test split as baseline for fair comparison
set.seed(101422)
biomarker_split_lasso <- biomarker_union %>%
    initial_split(prop = 0.8)

# Prepare matrices for glmnet
train_data <- training(biomarker_split_lasso)
test_data <- testing(biomarker_split_lasso)

x_train <- train_data %>%
    select(-class) %>%
    as.matrix()
y_train <- train_data %>%
    pull(class) %>%
    as.numeric()

x_test <- test_data %>%
    select(-class) %>%
    as.matrix()
y_test <- test_data %>%
    pull(class) %>%
    as.numeric()

# Cross-validation to find optimal lambda
cat("Performing cross-validation to select optimal lambda...\n")
set.seed(101422)
cv_lasso <- cv.glmnet(x_train, y_train,
    family = "binomial",
    alpha = 1, # alpha = 1 is LASSO
    nfolds = 10,
    type.measure = "class"
)

cat("Optimal lambda (min CV error):", cv_lasso$lambda.min, "\n")
cat("Optimal lambda (1SE rule):", cv_lasso$lambda.1se, "\n\n")

# Fit LASSO with optimal lambda (using 1SE rule for more regularization/simplicity)
lasso_fit <- glmnet(x_train, y_train,
    family = "binomial",
    alpha = 1,
    lambda = cv_lasso$lambda.1se
)

# Extract non-zero coefficients (selected proteins)
lasso_coefs <- coef(lasso_fit, s = cv_lasso$lambda.1se)
selected_proteins <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
selected_proteins <- selected_proteins[selected_proteins != "(Intercept)"]

cat("LASSO Selected Proteins:\n")
cat("Number of proteins:", length(selected_proteins), "\n")
print(selected_proteins)
cat("\n")

cat("LASSO Coefficients:\n")
lasso_coef_df <- data.frame(
    protein = rownames(lasso_coefs)[which(lasso_coefs != 0)],
    coefficient = as.numeric(lasso_coefs[which(lasso_coefs != 0)])
)
print(lasso_coef_df)
cat("\n")

# ===============================================================================
# STEP 4: EVALUATE LASSO PANEL PERFORMANCE
# ===============================================================================

cat("STEP 4: Evaluating LASSO Panel Performance\n")
cat(strrep("-", 70), "\n")

# Get predictions on test set
lasso_pred_prob <- predict(lasso_fit,
    newx = x_test,
    s = cv_lasso$lambda.1se,
    type = "response"
)[, 1]

lasso_pred_class <- ifelse(lasso_pred_prob > 0.5, 1, 0)

# Calculate metrics
lasso_results <- tibble(
    truth_class = factor(y_test, levels = c(0, 1)),
    pred = lasso_pred_prob,
    pred_class = factor(lasso_pred_class, levels = c(0, 1))
) %>%
    class_metrics(
        estimate = pred_class,
        truth = truth_class,
        pred,
        event_level = "second"
    )

cat("LASSO Panel Performance Metrics:\n")
print(lasso_results)
cat("\n\n")

# ===============================================================================
# STEP 5: TRY EVEN SIMPLER - FIT REGULAR LOGISTIC REGRESSION ON LASSO SELECTED
# ===============================================================================

cat("STEP 5: Refitting Standard Logistic Regression on LASSO-Selected Proteins\n")
cat(strrep("-", 70), "\n")

# Create dataset with only LASSO-selected proteins
biomarker_lasso_selected <- biomarker_clean %>%
    select(group, any_of(selected_proteins)) %>%
    mutate(class = (group == "ASD")) %>%
    select(-group)

# Use same split
set.seed(101422)
biomarker_split_final <- biomarker_lasso_selected %>%
    initial_split(prop = 0.8)

# Fit standard logistic regression
fit_lasso_glm <- glm(class ~ .,
    data = training(biomarker_split_final),
    family = "binomial"
)

cat("Standard GLM Coefficients (on LASSO-selected proteins):\n")
print(summary(fit_lasso_glm)$coefficients)
cat("\n")

# Evaluate
lasso_glm_results <- testing(biomarker_split_final) %>%
    add_predictions(fit_lasso_glm, type = "response") %>%
    mutate(
        pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE)),
        truth_class = factor(class, levels = c(FALSE, TRUE))
    ) %>%
    class_metrics(
        estimate = pred_class,
        truth = truth_class,
        pred,
        event_level = "second"
    )

cat("Standard GLM Performance (on LASSO-selected proteins):\n")
print(lasso_glm_results)
cat("\n\n")

# ===============================================================================
# STEP 6: COMPARISON SUMMARY
# ===============================================================================

cat(strrep("=", 70), "\n")
cat("FINAL COMPARISON SUMMARY\n")
cat(strrep("=", 70), "\n\n")

comparison_summary <- bind_rows(
    baseline_results %>% mutate(
        method = "Baseline (Intersection)",
        n_proteins = length(proteins_baseline)
    ),
    lasso_results %>% mutate(
        method = "LASSO Panel",
        n_proteins = length(selected_proteins)
    ),
    lasso_glm_results %>% mutate(
        method = "Standard GLM on LASSO",
        n_proteins = length(selected_proteins)
    )
) %>%
    select(method, n_proteins, .metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)

cat("Performance Comparison Table:\n")
print(comparison_summary)
cat("\n")

# Calculate reduction in complexity
cat("Complexity Reduction:\n")
cat("Baseline panel size:", length(proteins_baseline), "proteins\n")
cat("LASSO panel size:", length(selected_proteins), "proteins\n")
cat(
    "Reduction:", length(proteins_baseline) - length(selected_proteins),
    sprintf(
        "(%.1f%% fewer)",
        100 * (1 - length(selected_proteins) / length(proteins_baseline))
    ), "\n"
)
cat("\n")

# Save results
save(proteins_baseline,
    proteins_union,
    selected_proteins,
    lasso_fit,
    cv_lasso,
    comparison_summary,
    file = "results/q4_simpler_panel_results.RData"
)

cat("Results saved to: results/q4_simpler_panel_results.RData\n")
cat("\nAnalysis complete!\n")
