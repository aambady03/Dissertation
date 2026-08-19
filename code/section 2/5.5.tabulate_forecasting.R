library(dplyr)
library(knitr)
library(ENMeval)
library(terra)

base_dir <- "/home/anaga-ambady/Documents/DISSERT/DATA/part 2. india"

# =========================================================================
# 1. TABLE 1: ENMeval Tuning Results Across All Candidates
# =========================================================================
results_raw <- eval.results(e)

tuning_table <- results_raw %>%
  select(
    `FC`                       = fc,
    `RM`                       = rm,
    `AICc`                     = AICc,
    `Delta AICc`               = delta.AICc,
    `Parameters (K)`          = ncoef,
    `Validation AUC`          = auc.val.avg,
    `10% Omission Rate`       = or.10p.avg,
    `CBI`                      = cbi.val.avg
  ) %>%
  arrange(`Delta AICc`) %>%
  mutate(
    `AICc`               = round(`AICc`, 2),
    `Delta AICc`         = round(`Delta AICc`, 2),
    `Validation AUC`    = round(`Validation AUC`, 3),
    `10% Omission Rate` = round(`10% Omission Rate`, 3),
    `CBI`               = round(`CBI`, 3)
  )

cat("=== TABLE 1: ENMeval Tuning Results ===\n")
print(kable(tuning_table, format = "simple"))

write.csv(
  tuning_table,
  file.path(base_dir, "table1_tuning_results.csv"),
  row.names = FALSE
)

# =========================================================================
# AUTOMATED MODEL SELECTION (Multi-Tiered Validation Rule)
# =========================================================================
# Avoids overfitted models (e.g., Delta AICc = 0 with low AUC & high Omission Rate)
# Strategy: 1) Omission Rate <= 0.20, 2) Validation AUC > 0.55, 3) Min Delta AICc
valid_candidates <- results_raw %>%
  filter(or.10p.avg <= 0.20, auc.val.avg > 0.55) %>%
  arrange(delta.AICc)

if (nrow(valid_candidates) > 0) {
  best_row <- valid_candidates[1, ]
  cat("\n[Selection Rule]: Optimal model selected meeting Validation AUC > 0.55 & 10% Omission Rate <= 0.20.\n")
} else {
  warning("No models met strict OR <= 0.20 & AUC > 0.55. Selecting candidate with highest Continuous Boyce Index (CBI).")
  best_row <- results_raw %>% arrange(desc(cbi.val.avg)) %>% slice(1)
}

# =========================================================================
# 2. FINAL MODEL SPECIFICATIONS & ACTIVE COEFFICIENTS
# =========================================================================
cat("\n=== FINAL MODEL SPECIFICATIONS ===\n")
cat("Selected Model Code   :", best_row$tune.args, "\n")
cat("Feature Classes       :", best_row$fc, "\n")
cat("Regularization (RM)   :", best_row$rm, "\n")
cat("Delta AICc            :", round(best_row$delta.AICc, 2), "\n")
cat("Validation AUC        :", round(best_row$auc.val.avg, 3), "\n")
cat("10% Omission Rate     :", round(best_row$or.10p.avg, 3), "\n")
cat("Continuous Boyce (CBI):", round(best_row$cbi.val.avg, 3), "\n\n")

best_mod <- readRDS(file.path(base_dir, "best_maxnet_model.rds"))
cat("=== ACTIVE FEATURE COEFFICIENTS (K =", length(best_mod$betas), ") ===\n")
print(round(best_mod$betas, 4))

# =========================================================================
# 3. VARIABLE IMPORTANCE SUMMARY
# =========================================================================
var_imp <- eval.variable.importance(e)[[best_row$tune.args]]

if (!is.null(var_imp)) {
  var_imp_table <- var_imp %>%
    rename(
      `Variable`               = variable,
      `Percent Contribution`   = percent.contribution,
      `Permutation Importance` = permutation.importance
    ) %>%
    mutate(
      `Percent Contribution`   = round(`Percent Contribution`, 2),
      `Permutation Importance` = round(`Permutation Importance`, 2)
    ) %>%
    arrange(desc(`Permutation Importance`))
} else {
  cat("\nCalculating manual Permutation Importance for maxnet model...\n")
  eval_coords <- rbind(occs[, c("lon", "lat")], bg_points[, c("lon", "lat")])
  eval_data   <- na.omit(terra::extract(envs, eval_coords)[, names(envs)])
  base_pred   <- predict(best_mod, eval_data, type = "cloglog")
  
  var_names <- names(envs)
  imp_list  <- numeric(length(var_names))
  names(imp_list) <- var_names
  
  set.seed(42)
  for (v in var_names) {
    perm_data <- eval_data
    perm_data[[v]] <- sample(perm_data[[v]])
    perm_pred <- predict(best_mod, perm_data, type = "cloglog")
    imp_list[v] <- sqrt(mean((base_pred - perm_pred)^2))
  }
  
  imp_pct <- round((imp_list / sum(imp_list)) * 100, 2)
  var_imp_table <- data.frame(
    `Variable`               = names(imp_pct),
    `Permutation Importance` = imp_pct,
    check.names              = FALSE
  )
}

cat("\n=== VARIABLE IMPORTANCE SUMMARY ===\n")
print(kable(var_imp_table, format = "simple"))

write.csv(
  var_imp_table,
  file.path(base_dir, "final_model_var_importance.csv"),
  row.names = FALSE
)