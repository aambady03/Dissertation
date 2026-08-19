# ==============================================================================
# SCRIPT 03: MODEL COMPARISON & SELECTION TABLE GENERATION
# ==============================================================================

library(dplyr)
library(purrr)
library(mgcv)

model_dir <- "~/Documents/DISSERT/EDITS/section 1/gam_rds"
required_models <- c("gam_gdd_k_200", "gam_gdd_k_100", "gam_gdd_xalt",
                     "gam_gdd_fs", "gam_gdd_noroad", "gam_gdd_noroad_sm", "gam_gdd_no_loc")

for (mod_name in required_models) {
  if (!exists(mod_name, envir = .GlobalEnv)) {
    assign(mod_name, readRDS(file.path(model_dir, paste0(mod_name, ".rds"))), envir = .GlobalEnv)
  }
}

get_spatial_k_check <- function(model) {
  kc <- tryCatch(k.check(model), error = function(e) NULL)
  if (is.null(kc)) return("—")
  sp_row <- rownames(kc)[grep("longitude", rownames(kc))]
  if (length(sp_row) == 0) return("—")
  
  k_idx <- round(kc[sp_row, "k-index"], 2)
  p_val <- kc[sp_row, "p-value"]
  p_str <- if (p_val < 0.001) "< 0.001" else sprintf("%.3f", p_val)
  return(paste0(k_idx, " (p = ", p_str, ")"))
}

# --- Panel A: Spatial & Interaction Structure (REML) ---
reml_models <- list("M1_noloc" = gam_gdd_no_loc, "M2_fs" = gam_gdd_fs, 
                    "M3_xalt" = gam_gdd_xalt, "M4_k100" = gam_gdd_k_100, "M5_k200" = gam_gdd_k_200)
reml_ref    <- list("M1_noloc" = "M5_k200", "M2_fs" = "M4_k100", 
                    "M3_xalt" = "M4_k100", "M4_k100" = "M5_k200", "M5_k200" = NA)

panel_a <- imap_dfr(reml_models, function(m, id) {
  ref_id <- reml_ref[[id]]
  if (is.na(ref_id)) {
    p_str <- "Ref"; delta_aic <- 0; ref_label <- "—"
  } else {
    ref_model <- reml_models[[ref_id]]
    lrt   <- tryCatch(anova(m, ref_model, test = "Chisq"), error = function(e) NULL)
    p_lrt <- if (!is.null(lrt) && nrow(lrt) > 1) lrt$`Pr(>Chi)`[2] else NA
    p_str <- ifelse(is.na(p_lrt), "N/A", ifelse(p_lrt < 0.001, "< 0.001", sprintf("%.3f", p_lrt)))
    delta_aic <- round(AIC(m) - AIC(ref_model), 2)
    ref_label <- ref_id
  }
  data.frame(
    Panel = "Panel A: REML", Model_ID = id, Reference = ref_label,
    Total_EDF = round(sum(m$edf), 2), Deviance_Expl_Pct = round(summary(m)$dev.expl * 100, 2),
    AIC = round(AIC(m), 2), delta_AIC = delta_aic, LRT_p_val = p_str, Spatial_k_check = get_spatial_k_check(m)
  )
})

# --- Panel B: Fixed Effects Sensitivity (ML) ---
ml_models <- list(
  "M6_noroad_sm" = update(gam_gdd_noroad_sm, method = "ML"),
  "M7_noroad"    = update(gam_gdd_noroad,    method = "ML"),
  "M8_full_ml"   = update(gam_gdd_k_100,     method = "ML")
)
ml_ref <- list("M6_noroad_sm" = "M7_noroad", "M7_noroad" = "M8_full_ml", "M8_full_ml" = NA)

panel_b <- imap_dfr(ml_models, function(m, id) {
  ref_id <- ml_ref[[id]]
  if (is.na(ref_id)) {
    p_str <- "Ref"; delta_aic <- 0; ref_label <- "—"
  } else {
    ref_model <- ml_models[[ref_id]]
    lrt   <- tryCatch(anova(m, ref_model, test = "Chisq"), error = function(e) NULL)
    p_lrt <- if (!is.null(lrt) && nrow(lrt) > 1) lrt$`Pr(>Chi)`[2] else NA
    p_str <- ifelse(is.na(p_lrt), "N/A", ifelse(p_lrt < 0.001, "< 0.001", sprintf("%.3f", p_lrt)))
    delta_aic <- round(AIC(m) - AIC(ref_model), 2)
    ref_label <- ref_id
  }
  data.frame(
    Panel = "Panel B: ML", Model_ID = id, Reference = ref_label,
    Total_EDF = round(sum(m$edf), 2), Deviance_Expl_Pct = round(summary(m)$dev.expl * 100, 2),
    AIC = round(AIC(m), 2), delta_AIC = delta_aic, LRT_p_val = p_str, Spatial_k_check = get_spatial_k_check(m)
  )
})

# Export Combined Table
dissertation_table <- bind_rows(panel_a, panel_b)
tables_dir <- "~/Documents/DISSERT/EDITS/section 1/results/tables"
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(dissertation_table, file.path(tables_dir, "final_model_selection_table.csv"), row.names = FALSE)