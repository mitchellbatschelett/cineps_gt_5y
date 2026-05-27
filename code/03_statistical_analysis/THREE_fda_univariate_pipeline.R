################################################################################
#
#   FDA UNIVARIATE EXPOSURE-METRIC PIPELINE
#   Fully-adjusted PFFR models of graph theory metric(density) curves on a
#   single neonatal exposure plus forced covariates. Run within the VPT group
#   to motivate stability selection (Methods 2.7.4, Supp Note 5).
#
#   PURPOSE
#     Fits one (metric, exposure) cell of the 28-cell univariate grid:
#     7 candidate neonatal exposures x 4 primary graph theory metrics.
#     All 28 cells are fully adjusted with the same forced covariate set.
#     Provides per-predictor LR tests and 1000-iteration bootstrap CIs on
#     the functional beta coefficient for the exposure of interest.
#
#   INPUTS
#     - data/analysis_ready/cohort_171VPT_postVQC.xlsx
#         (171 VPT only, clinical + GT metrics)
#
#   OUTPUTS  (per run, written to {output_dir})
#     - {metric}_univariate_{exposure}_FDA_results.RData
#     - {metric}_univariate_{exposure}_summary.csv
#     - {metric}_univariate_{exposure}_predictor_tests.csv
#     - {metric}_univariate_{exposure}_pffr_diagnostics.png
#     - {metric}_univariate_{exposure}_pffr_residuals.png
#
#   REQUIRES
#     R 4.4.0
#     refund, mgcv, fda, fda.usc, readxl, tidyverse, ggplot2, viridis,
#     gridExtra, parallel
#
#   COMPUTATIONAL COST
#     On a 24-core HPC node with 128 GB RAM, expect
#     ~1-6 hours per cell including the 1000-iteration bootstrap.
#
#   NOTES
#     - Methods 2.7.4: "Univariate PFFR models for SW were fit under the
#       gaussian family for computational tractability." All 28 cells in
#       this grid use gaussian(), including the SW row.
#     - Forced covariates (Methods 2.7.4): eTIV, sex, sriskscore,
#       age_at_5y_mri, Rel_Motion.
#     - Standardization: ALL continuous predictors are z-scored before
#       fitting -- both the continuous exposures of interest (bw_z, ga,
#       globalbrainscore2 post-sqrt, dwma_percent) and the continuous
#       forced covariates (eTIV, sriskscore, age_at_5y_mri, Rel_Motion).
#       Binary variables (sex, bpd2, anyrop, sepsis2) are NOT standardized.
#       This is encoded via vars_to_skip_scaling = binary_vars.
#     - GBA (globalbrainscore2) is sqrt-transformed in code (Methods 2.7.4).
#     - Variance explained reported as R^2 (gaussian; i.e.
#       summary()$r.sq).
#
#   USAGE
#     Driven by THREE_fda_univariate_runs.R, which sets BATCH_* variables
#     and sources this script. Standalone use: edit the CONFIGURATION
#     block below.
#
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

if (exists("BATCH_MODE") && BATCH_MODE == TRUE) {
  
  metric                <- BATCH_metric
  exposure              <- BATCH_exposure
  forced_covariates     <- BATCH_forced_covariates
  vars_to_skip_scaling  <- BATCH_vars_to_skip_scaling
  density_min           <- BATCH_density_min
  density_max           <- BATCH_density_max
  n_bootstrap           <- BATCH_n_bootstrap
  p_threshold           <- BATCH_p_threshold
  pffr_k_basis          <- BATCH_pffr_k_basis
  pffr_family           <- BATCH_pffr_family
  subjects_to_exclude   <- BATCH_subjects_to_exclude
  apply_gba_sqrt        <- BATCH_apply_gba_sqrt
  output_dir            <- BATCH_output_dir
  data_path             <- BATCH_data_path
  n_cores               <- BATCH_n_cores
  
  cat("\n*** RUNNING IN BATCH MODE ***\n\n")
  
} else {
  
  # ----------------------------------------------------------------------------
  # STANDALONE CONFIGURATION - edit for interactive single-run use
  # ----------------------------------------------------------------------------
  
  metric    <- ""
  exposure  <- ""   # one of the 7 candidates
  forced_covariates <- c("eTIV", "sex", "sriskscore",
                         "age_at_5y_mri", "Rel_Motion")
  
  # Skip scaling for binary variables only. All continuous predictors
  # (exposures of interest + continuous forced covariates) get z-scored.
  vars_to_skip_scaling <- c("sex", "bpd2", "anyrop", "sepsis2")
  
  density_min  <- 11
  density_max  <- 100
  pffr_k_basis <- 20
  library(mgcv)
  pffr_family  <- gaussian()
  n_bootstrap  <- 1000
  # --- Cores: pick ONE line below (comment out the other) ---
  # n_cores <- 24                                 # fixed count (e.g. cluster node)
  n_cores <- max(1, parallel::detectCores() - 1)  # auto-detect (laptop-safe default)
  p_threshold  <- 0.001
  subjects_to_exclude <- c()   # set per metric
  apply_gba_sqrt      <- TRUE  # FALSE for non-GBA exposures
  
  # Repo root located automatically (standalone mode); a fresh clone runs without edits.
  # Manual override if needed: set `repo_root <- "/full/path/to/repository"` before source().
  .is_repo_root <- function(p) {
    dir.exists(file.path(p, "data", "analysis_ready")) && dir.exists(file.path(p, "code"))
  }
  .find_repo_root_from_path <- function(p) {
    p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    for (i in 1:8) {
      if (.is_repo_root(p)) return(p)
      parent <- dirname(p); if (parent == p) break; p <- parent
    }
    NULL
  }
  if (exists("repo_root", inherits = TRUE) && is.character(repo_root) &&
      length(repo_root) == 1 &&
      .is_repo_root(normalizePath(repo_root, winslash = "/", mustWork = FALSE))) {
    repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  } else {
    repo_root <- .find_repo_root_from_path(getwd())
    if (is.null(repo_root))
      stop("Could not locate repo root. Set it manually before source()-ing:\n",
           "  repo_root <- \"/full/path/to/repository\"")
  }
  cat(sprintf("Repo root: %s\n", repo_root))

  data_path       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
  output_dir_base <- file.path(repo_root, "results/fda_univariate")
  
  # Output: results/fda_univariate/{metric}/{metric}_univariate_{exposure}_11-100/
  density_str <- sprintf("%.0f-%.0f", density_min, density_max)
  output_dirname <- paste(metric, "univariate", exposure, density_str, sep = "_")
  output_dir <- file.path(output_dir_base, metric, output_dirname)
}

# ==============================================================================
# SETUP
# ==============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

required_packages <- c("refund", "fda", "mgcv", "readxl", "tidyverse",
                       "ggplot2", "viridis", "gridExtra", "fda.usc", "parallel")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages, repos = "https://cloud.r-project.org/",
                   lib = Sys.getenv("R_LIBS_USER"))
}

suppressPackageStartupMessages({
  library(refund);    library(fda);     library(mgcv);     library(readxl)
  library(tidyverse); library(ggplot2); library(viridis);  library(gridExtra)
  library(fda.usc);   library(parallel)
})

set.seed(42)

use_F_test <- inherits(pffr_family, "family") && pffr_family$family == "gaussian"

# All predictors in the model: exposure first, then forced covariates
all_predictors <- c(exposure, forced_covariates)

cat("--- Configuration Summary ---\n")
cat(sprintf("  Metric:               %s\n", metric))
cat(sprintf("  Exposure:             %s\n", exposure))
cat(sprintf("  Forced covariates:    %s\n",
            paste(forced_covariates, collapse = ", ")))
cat(sprintf("  Skip scaling for:     %s\n",
            paste(vars_to_skip_scaling, collapse = ", ")))
cat(sprintf("  Density range:        %.0f%% to %.0f%%\n",
            density_min, density_max))
cat(sprintf("  Family:               %s\n", pffr_family$family))
cat(sprintf("  Apply GBA sqrt:       %s\n", apply_gba_sqrt))
cat(sprintf("  Output:               %s\n", output_dir))

# ==============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 1: DATA LOADING AND PREPARATION\n")
cat("==========================================\n\n")

df <- read_excel(data_path)
n_original <- nrow(df)
cat(sprintf("Loaded data: %d subjects x %d variables\n", nrow(df), ncol(df)))

# --- 1.1 Metric-specific outlier exclusion (Methods 2.7.2) ---
if (length(subjects_to_exclude) > 0) {
  n_before <- nrow(df)
  df <- df[!df$ID %in% subjects_to_exclude, ]
  cat(sprintf("Excluded %d metric-specific outlier(s) by ID (%s). N = %d\n",
              n_before - nrow(df),
              paste(subjects_to_exclude, collapse = ", "), nrow(df)))
}

# --- 1.2 GBA sqrt transform (Methods 2.7.4) ---
if (apply_gba_sqrt && exposure == "globalbrainscore2") {
  cat(sprintf("\nApplying sqrt transform to globalbrainscore2\n"))
  cat(sprintf("  Before: range [%.1f, %.1f], median %.1f\n",
              min(df$globalbrainscore2, na.rm = TRUE),
              max(df$globalbrainscore2, na.rm = TRUE),
              median(df$globalbrainscore2, na.rm = TRUE)))
  df$globalbrainscore2 <- sqrt(df$globalbrainscore2)
  cat(sprintf("  After:  range [%.2f, %.2f], median %.2f\n",
              min(df$globalbrainscore2, na.rm = TRUE),
              max(df$globalbrainscore2, na.rm = TRUE),
              median(df$globalbrainscore2, na.rm = TRUE)))
}

# --- 1.3 Extract functional response ---
metric_cols <- grep(paste0("^", metric, "_"), colnames(df), value = TRUE)
cat(sprintf("\nFound %d density columns for %s\n", length(metric_cols), metric))
if (length(metric_cols) == 0) {
  stop(paste0("No columns matching '", metric, "_*'."))
}

densities_all  <- as.numeric(gsub(paste0(metric, "_"), "", metric_cols))
keep_densities <- densities_all >= density_min & densities_all <= density_max
densities      <- densities_all[keep_densities]
metric_cols    <- metric_cols[keep_densities]

cat(sprintf("Density range: %.2f%% to %.2f%% (%d points)\n",
            min(densities), max(densities), length(densities)))

Y <- as.matrix(df[, metric_cols])
rownames(Y) <- df$ID
colnames(Y) <- densities

# --- 1.4 Prepare predictors ---
predictor_df <- df %>%
  dplyr::select(ID, all_of(all_predictors)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(.)))

cat("\nMissingness in predictors:\n")
any_missing <- FALSE
for (var in all_predictors) {
  n_miss <- sum(is.na(predictor_df[[var]]))
  if (n_miss > 0) {
    cat(sprintf("  %s: %d missing (%.1f%%)\n",
                var, n_miss, 100 * n_miss / nrow(predictor_df)))
    any_missing <- TRUE
  }
}
if (!any_missing) cat("  No missing values\n")

complete_cases        <- complete.cases(predictor_df)
Y_complete            <- Y[complete_cases, ]
predictor_df_complete <- predictor_df[complete_cases, ]

cat(sprintf("\nComplete cases: %d / %d\n",
            sum(complete_cases), nrow(predictor_df)))

# --- 1.5 Standardize continuous variables ---
continuous_vars <- setdiff(all_predictors, vars_to_skip_scaling)

predictor_df_scaled <- predictor_df_complete
for (var in continuous_vars) {
  if (var %in% colnames(predictor_df_scaled)) {
    predictor_df_scaled[[var]] <- scale(predictor_df_scaled[[var]])[, 1]
  }
}

cat("\nContinuous predictors standardized (mean = 0, SD = 1):\n")
cat("  ", if (length(continuous_vars) > 0)
  paste(continuous_vars, collapse = ", ") else "(none)", "\n")
cat("Predictors not standardized:\n")
cat("  ", paste(vars_to_skip_scaling, collapse = ", "), "\n")

# Aliases used downstream
Y_clean                   <- Y_complete
predictor_df_clean        <- predictor_df_complete
predictor_df_scaled_clean <- predictor_df_scaled

cat(sprintf("\nFinal analysis sample: %d subjects\n", nrow(Y_clean)))

# ==============================================================================
# SECTION 2: PFFR MODEL FITTING
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 2: PFFR MODEL FITTING\n")
cat("==========================================\n\n")

pffr_data   <- predictor_df_scaled_clean
pffr_data$Y <- Y_clean

formula_parts <- paste(all_predictors, collapse = " + ")
pffr_formula  <- as.formula(paste("Y ~", formula_parts))

cat("Fitting pffr() model...\n")
cat(sprintf("  Formula: Y ~ %s\n", formula_parts))
cat(sprintf("  Family: %s\n", pffr_family$family))
cat("  Basis: P-splines with k = 20, first-order difference penalty\n\n")

pffr_fit <- pffr(
  pffr_formula,
  yind      = densities,
  data      = pffr_data,
  bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
  family    = pffr_family
)
cat("pffr() model fitted successfully.\n\n")
cat("Model Summary:\n==============\n")
print(summary(pffr_fit))

# ==============================================================================
# SECTION 3: HYPOTHESIS TESTING
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 3: HYPOTHESIS TESTING\n")
cat("==========================================\n\n")

# --- 3.1 Variance Explained ---
get_var_explained <- function(model) {
  s <- summary(model)
  if (use_F_test) s$r.sq else s$dev.expl
}

pffr_intercept_only <- pffr(
  Y ~ 1, yind = densities, data = pffr_data,
  bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
  family    = pffr_family
)

ve_full      <- get_var_explained(pffr_fit)
ve_intercept <- get_var_explained(pffr_intercept_only)
functional_ve <- (ve_full - ve_intercept) / (1 - ve_intercept)
ve_label     <- if (use_F_test) "R^2" else "Deviance explained"

cat(sprintf("--- Variance Explained (%s) ---\n\n", ve_label))
cat(sprintf("Full model:                   %.4f (%.2f%%)\n",
            ve_full, ve_full * 100))
cat(sprintf("Functional %-9s:        %.4f (%.2f%%)\n\n",
            ve_label, functional_ve, functional_ve * 100))

# --- 3.2 Overall LR Test ---
cat("--- Overall Likelihood Ratio Test ---\n\n")

constant_parts   <- paste0("c(", all_predictors, ")", collapse = " + ")
constant_formula <- as.formula(paste("Y ~ 1 +", constant_parts))

pffr_constant <- pffr(
  constant_formula, yind = densities, data = pffr_data,
  bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
  family    = pffr_family
)

mtest_full     <- pffr_fit
mtest_constant <- pffr_constant
class(mtest_full)     <- class(mtest_full)[-1]
class(mtest_constant) <- class(mtest_constant)[-1]

lr_test_overall <- anova(mtest_constant, mtest_full,
                         test = if (use_F_test) "F" else "Chisq")

if (use_F_test) {
  cat(sprintf("  F        = %.2f\n", lr_test_overall$F[2]))
  cat(sprintf("  p-value  = %.2e\n", lr_test_overall$`Pr(>F)`[2]))
  overall_p <- lr_test_overall$`Pr(>F)`[2]
} else {
  cat(sprintf("  Deviance = %.4f\n", lr_test_overall$Deviance[2]))
  cat(sprintf("  p-value  = %.2e\n", lr_test_overall$`Pr(>Chi)`[2]))
  overall_p <- lr_test_overall$`Pr(>Chi)`[2]
}
overall_model_significant <- overall_p < p_threshold
cat(sprintf("  Significant at p < %.3f: %s\n\n", p_threshold,
            ifelse(overall_model_significant, "YES", "NO")))

# --- 3.3 Per-Predictor LR Tests ---
cat("--- Per-Predictor Likelihood Ratio Tests ---\n\n")

predictor_tests <- data.frame(
  Predictor      = character(),
  Test_Statistic = numeric(),
  p_value        = numeric(),
  functional_ve  = numeric(),
  semi_partial_r = numeric(),
  stringsAsFactors = FALSE
)

for (pred in all_predictors) {
  other_preds <- setdiff(all_predictors, pred)
  reduced_formula <- as.formula(
    paste("Y ~", paste(other_preds, collapse = " + "))
  )
  
  pffr_reduced <- tryCatch(
    pffr(reduced_formula, yind = densities, data = pffr_data,
         bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
         family    = pffr_family),
    error = function(e) NULL
  )
  if (is.null(pffr_reduced)) {
    cat(sprintf("  %s: Reduced model fit failed; skipping.\n", pred)); next
  }
  mtest_reduced <- pffr_reduced
  class(mtest_reduced) <- class(mtest_reduced)[-1]
  
  lr_test <- tryCatch(
    anova(mtest_reduced, mtest_full,
          test = if (use_F_test) "F" else "Chisq"),
    error = function(e) NULL
  )
  if (is.null(lr_test)) {
    cat(sprintf("  %s: LR test failed; skipping.\n", pred)); next
  }
  
  single_formula <- as.formula(paste("Y ~ 1 +", pred))
  pffr_single <- tryCatch(
    pffr(single_formula, yind = densities, data = pffr_data,
         bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
         family    = pffr_family),
    error = function(e) NULL
  )
  func_ve_pred <- if (!is.null(pffr_single))
    (get_var_explained(pffr_single) - ve_intercept) / (1 - ve_intercept)
  else NA
  
  ve_reduced     <- get_var_explained(pffr_reduced)
  semi_partial_r <- sqrt(abs(ve_full - ve_reduced))
  
  pred_stat <- if (use_F_test) lr_test$F[2]        else lr_test$Deviance[2]
  pred_p    <- if (use_F_test) lr_test$`Pr(>F)`[2] else lr_test$`Pr(>Chi)`[2]
  
  predictor_tests <- rbind(predictor_tests, data.frame(
    Predictor      = pred,
    Test_Statistic = pred_stat,
    p_value        = pred_p,
    functional_ve  = func_ve_pred,
    semi_partial_r = semi_partial_r,
    stringsAsFactors = FALSE
  ))
  
  sig_marker <- ifelse(pred_p < p_threshold, " ***",
                       ifelse(pred_p < 0.01,        " **",
                              ifelse(pred_p < 0.05,        " *",  "")))
  cat(sprintf("  %-20s %s = %8.2f, p = %.4f, func_%s = %.4f, sr = %.4f%s\n",
              pred, if (use_F_test) "F" else "Dev",
              pred_stat, pred_p, ve_label, func_ve_pred,
              semi_partial_r, sig_marker))
}

cat(sprintf("\n  Significance: *** p < %.3f, ** p < 0.01, * p < 0.05\n",
            p_threshold))

write.csv(predictor_tests,
          file.path(output_dir,
                    paste0(metric, "_univariate_", exposure,
                           "_predictor_tests.csv")),
          row.names = FALSE)
cat(sprintf("\n  Saved: %s_univariate_%s_predictor_tests.csv\n",
            metric, exposure))

# ==============================================================================
# SECTION 4: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 4: BOOTSTRAP CONFIDENCE INTERVALS\n")
cat("==========================================\n\n")

cat(sprintf("Running %d bootstrap iterations on %d cores...\n",
            n_bootstrap, n_cores))

bootstrap_coefs <- tryCatch({
  coefboot.pffr(pffr_fit, B = n_bootstrap,
                ncpus = n_cores, parallel = "multicore")
}, error = function(e) {
  cat(sprintf("coefboot.pffr error: %s\n", e$message)); NULL
})

if (is.null(bootstrap_coefs)) {
  cat("Bootstrap CIs not available for this model.\n")
} else {
  cat("Bootstrap complete.\n")
}

# ==============================================================================
# SECTION 5: MODEL DIAGNOSTICS
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 5: MODEL DIAGNOSTICS\n")
cat("==========================================\n\n")

png(file.path(output_dir,
              paste0(metric, "_univariate_", exposure, "_pffr_diagnostics.png")),
    width = 12, height = 10, units = "in", res = 300)
pffr.check(pffr_fit)
dev.off()
cat(sprintf("  - %s_univariate_%s_pffr_diagnostics.png\n", metric, exposure))

Y_fitted       <- fitted(pffr_fit)
residuals_pffr <- Y_clean - Y_fitted

resid_df <- as.data.frame(residuals_pffr) %>%
  mutate(ID = predictor_df_clean$ID) %>%
  pivot_longer(cols = -ID, names_to = "density", values_to = "residual") %>%
  mutate(density = as.numeric(density))

p_resid <- ggplot(resid_df, aes(x = density, y = residual, group = ID)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Network Density (%)",
       y = "Residual (Observed - Fitted)",
       title = paste("Residual Curves -", metric, "~", exposure)) +
  theme_minimal()

ggsave(file.path(output_dir,
                 paste0(metric, "_univariate_", exposure, "_pffr_residuals.png")),
       p_resid, width = 10, height = 6, dpi = 300)
cat(sprintf("  - %s_univariate_%s_pffr_residuals.png\n", metric, exposure))

# ==============================================================================
# SECTION 6: SAVE RESULTS
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 6: SAVE RESULTS\n")
cat("==========================================\n\n")

summary_results <- data.frame(
  Metric                 = metric,
  Exposure               = exposure,
  Predictors             = paste(all_predictors, collapse = ", "),
  Family                 = pffr_family$family,
  Density_Range          = paste0(density_min, "-", density_max, "%"),
  N_Original             = n_original,
  N_Excluded_Outliers    = length(subjects_to_exclude),
  N_Subjects_Final       = nrow(Y_clean),
  Variance_Explained     = ve_full,
  Functional_VE          = functional_ve,
  VE_Label               = ve_label,
  Overall_Test_Statistic = if (use_F_test) lr_test_overall$F[2]
  else            lr_test_overall$Deviance[2],
  Overall_Test_Type      = if (use_F_test) "F" else "Chisq",
  Overall_p              = overall_p,
  Overall_Significant    = overall_model_significant,
  P_Threshold            = p_threshold,
  GBA_Sqrt_Applied       = (apply_gba_sqrt && exposure == "globalbrainscore2"),
  stringsAsFactors       = FALSE
)

write.csv(summary_results,
          file.path(output_dir,
                    paste0(metric, "_univariate_", exposure, "_summary.csv")),
          row.names = FALSE)
cat(sprintf("  - %s_univariate_%s_summary.csv\n", metric, exposure))

save(
  metric, exposure, densities, density_min, density_max, pffr_family,
  Y_clean, Y_complete, predictor_df_clean, predictor_df_scaled_clean,
  pffr_fit, pffr_intercept_only, pffr_constant,
  predictor_tests, summary_results,
  lr_test_overall, ve_full, ve_intercept, functional_ve, ve_label,
  bootstrap_coefs,
  subjects_to_exclude, all_predictors, vars_to_skip_scaling,
  apply_gba_sqrt,
  file = file.path(output_dir,
                   paste0(metric, "_univariate_", exposure, "_FDA_results.RData"))
)
cat(sprintf("  - %s_univariate_%s_FDA_results.RData\n", metric, exposure))

cat("\n\nSession Information:\n")
print(sessionInfo())

# ==============================================================================
# RUN SUMMARY
# ==============================================================================

cat("\n\n==========================================\n")
cat("RUN SUMMARY\n")
cat("==========================================\n\n")

cat(sprintf("Metric:       %s\n", metric))
cat(sprintf("Exposure:     %s\n", exposure))
cat(sprintf("Predictors:   %s\n", paste(all_predictors, collapse = ", ")))
cat(sprintf("Family:       %s\n", pffr_family$family))
cat(sprintf("N:            %d\n", nrow(Y_clean)))

cat(sprintf("\n%s (full):      %.4f\n", ve_label, ve_full))
cat(sprintf("Functional %s:  %.4f (%.2f%%)\n\n",
            ve_label, functional_ve, functional_ve * 100))

cat("Overall LR test:\n")
if (use_F_test) {
  cat(sprintf("  F = %.2f, p = %.2e\n", lr_test_overall$F[2], overall_p))
} else {
  cat(sprintf("  Chi^2 = %.4f, p = %.2e\n",
              lr_test_overall$Deviance[2], overall_p))
}
cat(sprintf("  Significant at p < %.3f: %s\n\n", p_threshold,
            ifelse(overall_model_significant, "YES", "NO")))

cat("Per-predictor tests:\n")
print(predictor_tests[order(predictor_tests$p_value), ], row.names = FALSE)

cat(sprintf("\nOutputs in:\n%s\n", output_dir))