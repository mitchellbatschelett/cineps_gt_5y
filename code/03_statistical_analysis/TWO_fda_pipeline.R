################################################################################
#
#   FDA PIPELINE: Function-on-Scalar Regression for VPT vs FT group differences
#   in global graph theory metrics across network density.
#
#   PURPOSE
#     Fits a PFFR model of metric(density) on Group plus optional covariates
#     for one (metric, covariate configuration) cell. Produces 1000-iteration
#     bootstrap confidence intervals on the functional beta coefficients.
#
#   INPUTS
#     - data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx
#         (clinical + GT metrics, 216 subjects, Group included)
#
#   OUTPUTS  (per run, written to {output_dir})
#     - {metric}_FDA_results.RData    Fitted pffr_fit, bootstrap_coefs, tests
#     - {metric}_summary.csv          One-row run summary (N, R^2/dev, p-values)
#     - {metric}_predictor_tests.csv  Per-predictor LR test results
#     - {metric}_pffr_diagnostics.png Q-Q, response vs fitted, residuals, k-check
#     - {metric}_pffr_residuals.png   Per-subject residual curves
#
#   REQUIRES
#     R 4.4.0
#     refund, mgcv, fda, fda.usc, readxl, tidyverse, ggplot2, viridis,
#     gridExtra, parallel
#
#   COMPUTATIONAL COST
#     These models are expensive. On a 48-core HPC node with 128 GB RAM,
#     gaussian-family runs (strength, normalized/raw GE, normalized/raw ACC)
#     take several hours each including the 1000-iteration bootstrap. SW
#     under the scat family can take over a week per configuration. The
#     submit_FDA_array_group_diff.lsf script submits all 36 runs in parallel
#     as an array job; total wall time is gated by the slowest scat-SW run.
#     Full reproduction from scratch is impractical on a workstation; the
#     saved *_FDA_results.RData files in results/fda_group_differences/ are
#     the recommended starting point for figure regeneration and inspection.
#
#   USAGE
#     Driven by TWO_fda_group_difference_runs.R, which sets BATCH_* variables
#     and sources this script. To run a single configuration interactively,
#     edit the STANDALONE CONFIGURATION block below and source this file.
#
#
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

if (exists("BATCH_MODE") && BATCH_MODE == TRUE) {
  
  # Driven by TWO_fda_group_difference_runs.R
  metric                 <- BATCH_metric
  all_predictors         <- BATCH_all_predictors
  vars_to_skip_scaling   <- BATCH_vars_to_skip_scaling
  density_min            <- BATCH_density_min
  density_max            <- BATCH_density_max
  n_bootstrap            <- BATCH_n_bootstrap
  p_threshold            <- BATCH_p_threshold
  pffr_k_basis           <- BATCH_pffr_k_basis
  pffr_family            <- BATCH_pffr_family
  subjects_to_exclude    <- BATCH_subjects_to_exclude
  output_dir             <- BATCH_output_dir
  data_path              <- BATCH_data_path
  n_cores                <- BATCH_n_cores
  
  cat("\n*** RUNNING IN BATCH MODE ***\n\n")
  
} else {
  
  # ----------------------------------------------------------------------------
  # STANDALONE CONFIGURATION - edit for interactive single-run use
  # ----------------------------------------------------------------------------
  
  # Metric column prefix (one of: str, rand_norm_wei_GE, rand_norm_wei_ACC,
  # rand_norm_wei_SW, GE, ACC). See data/intermediate/graph_theory_metrics_per_subject.xlsx.
  metric <- ""
  
  # Predictors (Group must be first for downstream figure scripts to find it)
  all_predictors <- c("Group", "age_at_5y_mri", "eTIV", "sex", "Rel_Motion")
  
  # Variables not to standardize. Group and sex are binary; everything else
  # in all_predictors is continuous and will be z-scored before fitting.
  vars_to_skip_scaling <- c("Group", "sex")
  
  # Density range (matches Methods 2.7.2)
  density_min <- 11
  density_max <- 100
  
  # PFFR settings (matches Methods 2.7.3)
  pffr_k_basis <- 20
  pffr_family  <- gaussian()      # use scat() for rand_norm_wei_SW
  
  # Bootstrap settings
  n_bootstrap <- 1000
  # --- Cores: pick ONE line below (comment out the other) ---
  # n_cores <- 24                                 # fixed count (e.g. cluster node)
  n_cores <- max(1, parallel::detectCores() - 1)  # auto-detect (laptop-safe default)
  
  # Inference threshold (Methods 2.7.3)
  p_threshold <- 0.001
  
  # Metric-specific outlier IDs (Methods 2.7.2). Pass an empty vector if none.
  # Lookups for reference:
  #   str:               c(128, 713)
  #   rand_norm_wei_GE:  c(659)
  #   rand_norm_wei_ACC: c()
  #   rand_norm_wei_SW:  c(309, 321, 8155)
  #   GE  (raw):         c(128, 713)
  #   ACC (raw):         c(128, 713)
  subjects_to_exclude <- c()
  
  # ----------------------------------------------------------------------------
  # PATHS - cluster defaults; uncomment the LOCAL OVERRIDE block for laptop use
  # ----------------------------------------------------------------------------
  
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

  data_path       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
  output_dir_base <- file.path(repo_root, "results/fda_group_differences")
  
  # Output directory name follows the convention
  # {metric}_FDA_{predictors_joined_by_underscore}_{dmin}-{dmax}_fullsample
  pred_str    <- paste(all_predictors, collapse = "_")
  density_str <- sprintf("%.0f-%.0f", density_min, density_max)
  output_dir  <- file.path(output_dir_base,
                           paste(metric, "FDA", pred_str, density_str,
                                 "fullsample", sep = "_"))
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
  library(refund);   library(fda);     library(mgcv);     library(readxl)
  library(tidyverse); library(ggplot2); library(viridis); library(gridExtra)
  library(fda.usc);  library(parallel)
})

set.seed(42)

# Family determines test statistic for anova.gam (F for gaussian, Chisq otherwise)
use_F_test <- inherits(pffr_family, "family") && pffr_family$family == "gaussian"

# ==============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 1: DATA LOADING AND PREPARATION\n")
cat("========================================\n\n")

cat(sprintf("Metric:               %s\n", metric))
cat(sprintf("Predictors:           %s\n", paste(all_predictors, collapse = ", ")))
cat(sprintf("Skip scaling for:     %s\n",
            if (length(vars_to_skip_scaling) > 0)
              paste(vars_to_skip_scaling, collapse = ", ") else "none"))
cat(sprintf("Density range:        %.2f%% to %.2f%%\n", density_min, density_max))
cat(sprintf("Family:               %s\n", pffr_family$family))
cat(sprintf("Inference threshold:  p < %.3f\n", p_threshold))
cat(sprintf("Output directory:     %s\n\n", output_dir))

df <- read_excel(data_path)
n_original <- nrow(df)

if (length(subjects_to_exclude) > 0) {
  n_before <- nrow(df)
  df <- df[!df$ID %in% subjects_to_exclude, ]
  cat(sprintf("Excluded %d metric-specific outlier(s) by ID (%s). N = %d\n",
              n_before - nrow(df),
              paste(subjects_to_exclude, collapse = ", "), nrow(df)))
}

cat(sprintf("Data loaded. Dimensions: %d subjects x %d variables\n",
            nrow(df), ncol(df)))

# ------------------------------------------------------------------------------
# 1.1 Extract Graph Theory Functional Data
# ------------------------------------------------------------------------------

metric_cols <- grep(paste0("^", metric, "_"), colnames(df), value = TRUE)
cat(sprintf("\nFound %d density threshold columns for %s\n",
            length(metric_cols), metric))

if (length(metric_cols) == 0) {
  all_density_cols <- colnames(df)[grep("_\\d+\\.\\d+$", colnames(df))]
  available_prefixes <- unique(gsub("_\\d+\\.\\d+$", "", all_density_cols))
  stop(paste0("No columns found matching pattern: ", metric, "_*\n",
              "Available metric prefixes: ",
              paste(available_prefixes, collapse = ", ")))
}

densities_all  <- as.numeric(gsub(paste0(metric, "_"), "", metric_cols))
keep_densities <- densities_all >= density_min & densities_all <= density_max
densities      <- densities_all[keep_densities]
metric_cols    <- metric_cols[keep_densities]

cat(sprintf("Trimmed density range: %.2f%% to %.2f%% (%d points)\n",
            min(densities), max(densities), length(densities)))

Y <- as.matrix(df[, metric_cols])
rownames(Y) <- df$ID
colnames(Y) <- densities

cat(sprintf("\nFunctional response matrix Y: %d subjects x %d density points\n",
            nrow(Y), ncol(Y)))
cat(sprintf("%s range: [%.4f, %.4f]\n", metric,
            min(Y, na.rm = TRUE), max(Y, na.rm = TRUE)))

# ------------------------------------------------------------------------------
# 1.2 Predictor Frame and Complete-Case Filter
# ------------------------------------------------------------------------------

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

complete_cases <- complete.cases(predictor_df)
Y_complete             <- Y[complete_cases, ]
predictor_df_complete  <- predictor_df[complete_cases, ]

cat(sprintf("\nComplete cases: %d / %d\n",
            sum(complete_cases), nrow(predictor_df)))

# ------------------------------------------------------------------------------
# 1.3 Standardize Continuous Predictors
# ------------------------------------------------------------------------------

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
cat("  ", if (length(vars_to_skip_scaling) > 0)
  paste(vars_to_skip_scaling, collapse = ", ") else "(none)", "\n")

# Aliases used downstream (preserved from upstream pipeline naming)
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

# ------------------------------------------------------------------------------
# 3.1 Variance Explained (R^2 for gaussian, deviance explained for scat)
# ------------------------------------------------------------------------------

# For non-gaussian families, summary()$r.sq is misleading (Wood et al. 2016).
# Use summary()$dev.expl ("deviance explained") instead. Methods 2.7.3 uses
# this convention for scat-family SW models.
get_var_explained <- function(model) {
  s <- summary(model)
  if (use_F_test) s$r.sq else s$dev.expl
}

pffr_intercept_only <- pffr(
  Y ~ 1, yind = densities, data = pffr_data,
  bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
  family    = pffr_family
)

ve_full         <- get_var_explained(pffr_fit)
ve_intercept    <- get_var_explained(pffr_intercept_only)
functional_ve   <- (ve_full - ve_intercept) / (1 - ve_intercept)
ve_label        <- if (use_F_test) "R^2" else "Deviance explained"

cat(sprintf("--- Variance Explained (%s) ---\n\n", ve_label))
cat(sprintf("Full model:           %.4f (%.2f%%)\n", ve_full, ve_full * 100))
cat(sprintf("Functional %-9s: %.4f (%.2f%%)\n\n",
            ve_label, functional_ve, functional_ve * 100))

# ------------------------------------------------------------------------------
# 3.2 Overall LR Test (Methods 2.7.3: varying vs. constant coefficients)
# ------------------------------------------------------------------------------

cat("--- Overall Likelihood Ratio Test ---\n\n")

constant_parts    <- paste0("c(", all_predictors, ")", collapse = " + ")
constant_formula  <- as.formula(paste("Y ~ 1 +", constant_parts))

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

# ------------------------------------------------------------------------------
# 3.3 Per-Predictor LR Tests
# ------------------------------------------------------------------------------

cat("--- Per-Predictor Likelihood Ratio Tests ---\n\n")

if (length(all_predictors) == 1) {
  
  cat("  Single-predictor model: per-predictor test is identical to overall test.\n\n")
  predictor_tests <- data.frame(
    Predictor      = all_predictors[1],
    Test_Statistic = if (use_F_test) lr_test_overall$F[2]
    else            lr_test_overall$Deviance[2],
    p_value        = overall_p,
    functional_ve  = functional_ve,
    semi_partial_r = sqrt(functional_ve),
    stringsAsFactors = FALSE
  )
  
} else {
  
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
    
    # Functional VE for this predictor in isolation
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
    
    # Semi-partial: sqrt of unique VE drop when this predictor is removed
    ve_reduced       <- get_var_explained(pffr_reduced)
    semi_partial_r   <- sqrt(abs(ve_full - ve_reduced))
    
    pred_stat <- if (use_F_test) lr_test$F[2]            else lr_test$Deviance[2]
    pred_p    <- if (use_F_test) lr_test$`Pr(>F)`[2]     else lr_test$`Pr(>Chi)`[2]
    
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
}

write.csv(predictor_tests,
          file.path(output_dir, paste0(metric, "_predictor_tests.csv")),
          row.names = FALSE)
cat(sprintf("\n  Saved: %s_predictor_tests.csv\n", metric))

# ==============================================================================
# SECTION 4: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 4: BOOTSTRAP CONFIDENCE INTERVALS\n")
cat("==========================================\n\n")

cat(sprintf("Running %d bootstrap iterations on %d cores...\n",
            n_bootstrap, n_cores))
cat("(Wall time scales with metric, family, and core count; see header.)\n\n")

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

# Combined diagnostic panel (Q-Q, response vs fitted, residuals, k-check histogram)
png(file.path(output_dir, paste0(metric, "_pffr_diagnostics.png")),
    width = 12, height = 10, units = "in", res = 300)
pffr.check(pffr_fit)
dev.off()
cat(sprintf("  - %s_pffr_diagnostics.png\n", metric))

# Per-subject residual curves
Y_fitted      <- fitted(pffr_fit)
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
       title = paste("Residual Curves -", metric)) +
  theme_minimal()

ggsave(file.path(output_dir, paste0(metric, "_pffr_residuals.png")),
       p_resid, width = 10, height = 6, dpi = 300)
cat(sprintf("  - %s_pffr_residuals.png\n", metric))

# ==============================================================================
# SECTION 6: SAVE RESULTS
# ==============================================================================

cat("\n\n==========================================\n")
cat("SECTION 6: SAVE RESULTS\n")
cat("==========================================\n\n")

summary_results <- data.frame(
  Metric                  = metric,
  Predictors              = paste(all_predictors, collapse = ", "),
  Family                  = pffr_family$family,
  Density_Range           = paste0(density_min, "-", density_max, "%"),
  N_Original              = n_original,
  N_Excluded_Outliers     = length(subjects_to_exclude),
  N_Subjects_Final        = nrow(Y_clean),
  Variance_Explained      = ve_full,
  Functional_VE           = functional_ve,
  VE_Label                = ve_label,
  Overall_Test_Statistic  = if (use_F_test) lr_test_overall$F[2]
  else            lr_test_overall$Deviance[2],
  Overall_Test_Type       = if (use_F_test) "F" else "Chisq",
  Overall_p               = overall_p,
  Overall_Significant     = overall_model_significant,
  P_Threshold             = p_threshold
)

write.csv(summary_results,
          file.path(output_dir, paste0(metric, "_summary.csv")),
          row.names = FALSE)
cat(sprintf("  - %s_summary.csv\n", metric))

save(
  metric, densities, density_min, density_max, pffr_family,
  Y_clean, Y_complete, predictor_df_clean, predictor_df_scaled_clean,
  pffr_fit, pffr_intercept_only, pffr_constant,
  predictor_tests, summary_results,
  lr_test_overall, ve_full, ve_intercept, functional_ve, ve_label,
  bootstrap_coefs,
  subjects_to_exclude, all_predictors, vars_to_skip_scaling,
  file = file.path(output_dir, paste0(metric, "_FDA_results.RData"))
)
cat(sprintf("  - %s_FDA_results.RData\n", metric))

cat("\n\nSession Information:\n")
print(sessionInfo())

# ==============================================================================
# RUN SUMMARY
# ==============================================================================

cat("\n\n==========================================\n")
cat("RUN SUMMARY\n")
cat("==========================================\n\n")

cat(sprintf("Metric:       %s\n", metric))
cat(sprintf("Predictors:   %s\n", paste(all_predictors, collapse = ", ")))
cat(sprintf("Family:       %s\n", pffr_family$family))
cat(sprintf("Density:      %.2f%% to %.2f%%\n", density_min, density_max))
cat(sprintf("N:            %d", nrow(Y_clean)))
if (length(subjects_to_exclude) > 0) {
  cat(sprintf(" (%d metric-specific outlier(s) excluded)\n",
              length(subjects_to_exclude)))
} else {
  cat(" (no exclusions)\n")
}

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