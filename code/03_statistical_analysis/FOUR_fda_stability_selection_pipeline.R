################################################################################
#
#   FDA STABILITY SELECTION PIPELINE
#   Group LASSO with stability selection for graph theory metric ~ neonatal
#   exposures, residualized on forced covariates.
#
#   PURPOSE
#     Identifies which of six candidate neonatal exposures (BPD, BWZ, GA, GBA,
#     ROP, sepsis) carry independent explanatory information for a graph
#     theory metric in the VPT group. For stably-selected exposures, fits a
#     full-sample PFFR refit with 1000-iteration bootstrap CIs on the
#     functional beta coefficients.
#
#   INPUTS
#     - data/analysis_ready/cohort_171VPT_postVQC.xlsx
#         (clinical + GT metrics, 171 VPT only)
#
#   OUTPUTS  (per run, written to {output_dir})
#     - {metric}_selection_matrix.{rds,csv}     100x6 binary selection matrix
#                                                (consumed by Supp Note 7)
#     - {metric}_stability_selection_summary.csv  Selection frequencies and
#                                                 stable/unstable flags (Table 3)
#     - {metric}_stability_selection_frequencies.png  Figures 3A / 4A / 5
#     - {metric}_predictor_tests.csv           Per-predictor LR tests in
#                                              full-sample refit
#     - {metric}_model_summary.csv             Run-level summary (R^2,
#                                              functional R^2, overall p)
#     - {metric}_final_beta_{predictor}.png    Bootstrap CI plots
#                                              (Figures 3B / 4B)
#     - {metric}_pffr_diagnostics.png          pffr.check output
#     - {metric}_final_pffr_residuals.png      Per-subject residual curves
#     - {metric}_stabsel_results.RData         Full workspace for downstream
#                                              figure / supplement scripts
#
#   REQUIRES
#     refund, fda, mgcv, readxl, tidyverse, ggplot2, viridis, gridExtra,
#     parallel, grpreg
#
#   COMPUTATIONAL COST
#     100 splits each running a
#     PFFR residualization + FPCA + cross-validated group LASSO. On a recent
#     laptop, expect roughly 60 minutes for stability selection per metric,
#     plus many hours for the full-sample 1000-iteration bootstrap if
#     any exposure is stably selected. Total wall time for all 6 published
#     runs (4 metrics + 2 ACC sensitivities) is on the order of days. 
#     The saved RData files in results/fda_stability_selection/ are the 
#     recommended starting point for figure regeneration.
#
#   USAGE
#     Driven by FOUR_fda_stability_selection_runs.R, which sets BATCH_*
#     variables and sources this script. To run a single configuration
#     interactively, edit the STANDALONE CONFIGURATION block below.
#
#   NOTES
#     - Methods 2.7.4: 100 splits, 60% selection / 40% holdout, 70% stability
#       threshold, 10-fold CV with lambda.min, FPCA PVE threshold 99.5%
#       (with min 3 / max 10 PCs), inverse-variance PC weighting,
#       square-root transform on GBA, alpha 0.001 in the full-sample
#       refit.
#     - Six candidate exposures: BPD (binary), BWZ (continuous), GA (continuous),
#       GBA (continuous, sqrt-transformed), ROP (binary), sepsis (binary).
#       DWMA was excluded prior to stability selection (Methods 2.7.4).
#     - Forced covariates (entered unpenalized into every model in B, D):
#       eTIV, sex, sriskscore, age_at_5y_mri, Rel_Motion.
#     - For SW: pffr_family_final = scat() per Methods 2.7.3/2.7.4 (heavy-
#       tailed residuals). All other metrics use gaussian().
#     - Sensitivity branches for ACC (Supp Note 6, Supp Figs 15-16):
#         "main"          : analysis as published in Figure 3
#         "high_gba_rem"  : exclude VPT participants with globalbrainscore2 > 14
#                           (yields the 4 participants reported in Supp Fig 15)
#         "gba_binary"    : replace globalbrainscore2 with globalcatmod
#                           (>=8 vs <8) and skip the sqrt transform
#       Branches are encoded by BATCH_sensitivity_branch.
#
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

if (exists("BATCH_MODE") && BATCH_MODE == TRUE) {
  
  # Driven by FOUR_fda_stability_selection_runs.R
  metric              <- BATCH_metric
  forced_covariates   <- BATCH_forced_covariates
  exposures           <- BATCH_exposures
  categorical_exposures <- BATCH_categorical_exposures
  categorical_forced  <- BATCH_categorical_forced
  density_min         <- BATCH_density_min
  density_max         <- BATCH_density_max
  pffr_family_final   <- BATCH_pffr_family_final
  subjects_to_exclude <- BATCH_subjects_to_exclude
  sensitivity_branch  <- BATCH_sensitivity_branch
  gbs2_transform      <- BATCH_gbs2_transform
  output_dir          <- BATCH_output_dir
  data_path           <- BATCH_data_path
  
  cat("\n*** RUNNING IN BATCH MODE ***\n\n")
  
} else {
  
  # ----------------------------------------------------------------------------
  # STANDALONE CONFIGURATION - edit for interactive single-run use
  # ----------------------------------------------------------------------------
  
  metric              <- ""
  forced_covariates   <- c("eTIV", "sex", "sriskscore", "age_at_5y_mri", "Rel_Motion")
  exposures           <- c("bpd2", "bw_z", "ga", "globalbrainscore2", "anyrop", "sepsis2")
  categorical_exposures <- c("bpd2", "anyrop", "sepsis2")
  categorical_forced    <- c("sex")
  
  density_min         <- 11
  density_max         <- 100
  
  # gaussian() for str, normalized/raw GE, normalized/raw ACC; scat() for SW
  library(mgcv)
  pffr_family_final   <- gaussian()
  
  # Metric-specific outlier IDs (Methods 2.7.2). Set per metric.
  subjects_to_exclude <- c()
  
  # Sensitivity branch: "main", "high_gba_rem", or "gba_binary"
  sensitivity_branch  <- ""
  
  # GBA transform: "sqrt" for primary analyses; "none" when sensitivity_branch
  # == "gba_binary" (since globalcatmod is already 0/1)
  gbs2_transform      <- "sqrt"
  
  # ----------------------------------------------------------------------------
  # PATHS
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

  data_path       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
  output_dir_base <- file.path(repo_root, "results/fda_stability_selection")
  
  # Output directory follows the existing convention:
  #   {metric}_stabsel_{dmin}-{dmax}/        for main analyses
  #   ACC_sensitivity/{metric}_stabsel_{dmin}-{dmax}_{branch}/  for sensitivities
  density_str <- sprintf("%.0f-%.0f", density_min, density_max)
  if (sensitivity_branch == "main") {
    output_dirname <- paste(metric, "stabsel", density_str, sep = "_")
    output_dir <- file.path(output_dir_base, output_dirname)
  } else {
    branch_suffix <- switch(sensitivity_branch,
                            "high_gba_rem" = "high_gba_rem",
                            "gba_binary"   = "gba_binary")
    output_dirname <- paste(metric, "stabsel", density_str,
                            branch_suffix, sep = "_")
    output_dir <- file.path(output_dir_base, "ACC_sensitivity", output_dirname)
  }
}

# ==============================================================================
# FIXED PARAMETERS (manuscript values; not run-to-run configurable)
# ==============================================================================

fpca_pve_threshold  <- 0.995    # Methods 2.7.4
fpca_min_npc        <- 3
fpca_max_npc        <- 10

n_splits            <- 100      # Methods 2.7.4
selection_frac      <- 0.60     # Methods 2.7.4
stability_threshold <- 0.70     # Methods 2.7.4
lasso_criterion     <- "lambda.min"  # Methods 2.7.4

pffr_k_basis        <- 20       # Methods 2.7.3 / 2.7.4
n_cv_folds          <- 10       # Methods 2.7.4

n_bootstrap         <- 1000     # Methods 2.7.4
p_threshold         <- 0.001    # Methods 2.7.4

master_seed         <- 42

# ==============================================================================
# SETUP
# ==============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

required_packages <- c("refund", "fda", "mgcv", "readxl", "tidyverse",
                       "ggplot2", "viridis", "gridExtra", "parallel", "grpreg")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages, repos = "https://cloud.r-project.org/")
}

suppressPackageStartupMessages({
  library(refund);    library(fda);     library(mgcv);     library(readxl)
  library(tidyverse); library(ggplot2); library(viridis);  library(gridExtra)
  library(parallel);  library(grpreg)
})

set.seed(master_seed)

# LR test type depends on family of the final inference model
lr_test_type <- if (inherits(pffr_family_final, "family") &&
                    pffr_family_final$family == "gaussian") "F" else "Chisq"

cat("--- Configuration Summary ---\n")
cat(sprintf("  Metric:              %s\n", metric))
cat(sprintf("  Density range:       %.0f%% to %.0f%%\n", density_min, density_max))
cat(sprintf("  Forced covariates:   %s\n", paste(forced_covariates, collapse = ", ")))
cat(sprintf("  Exposures:           %s\n", paste(exposures, collapse = ", ")))
cat(sprintf("  Sensitivity branch:  %s\n", sensitivity_branch))
cat(sprintf("  GBA transform:       %s\n", gbs2_transform))
cat(sprintf("  Stability selection: %d splits, %.0f%%/%.0f%%, threshold %.0f%%\n",
            n_splits, selection_frac * 100, (1 - selection_frac) * 100,
            stability_threshold * 100))
cat(sprintf("  Final model family:  %s\n",
            if (inherits(pffr_family_final, "general.family"))
              "scat" else pffr_family_final$family))
cat(sprintf("  LR test type:        %s\n", lr_test_type))
cat(sprintf("  Output:              %s\n", output_dir))

# ==============================================================================
# SECTION A: DATA LOADING AND PREPARATION
# ==============================================================================

cat("\n\n========================================\n")
cat("SECTION A: DATA LOADING AND PREPARATION\n")
cat("========================================\n\n")

df <- read_excel(data_path)
n_original <- nrow(df)
cat(sprintf("Loaded data: %d subjects x %d variables\n", nrow(df), ncol(df)))

# --- A.1 Sensitivity branch transforms (applied BEFORE complete-case filter) ---

if (sensitivity_branch == "high_gba_rem") {
  # Supp Note 6 / Supp Fig 15 sensitivity. Filters by globalbrainscore2 > 14;
  # in this cohort that yields the 4 participants reported in the manuscript.
  n_before <- nrow(df)
  df <- df[df$globalbrainscore2 <= 14 | is.na(df$globalbrainscore2), ]
  n_dropped <- n_before - nrow(df)
  cat(sprintf("Sensitivity 'high_gba_rem': removed %d participant(s) with globalbrainscore2 > 14. N = %d\n",
              n_dropped, nrow(df)))
} else if (sensitivity_branch == "gba_binary") {
  # Supp Note 6 / Supp Fig 16 sensitivity. Replace continuous globalbrainscore2
  # with the binary globalcatmod (1 = score >= 8, 0 = score < 8). The runs
  # driver should set gbs2_transform = "none" for this branch since the variable
  # is already binary.
  stopifnot("globalcatmod" %in% colnames(df))
  exposures[exposures == "globalbrainscore2"] <- "globalcatmod"
  if (!("globalcatmod" %in% categorical_exposures)) {
    categorical_exposures <- c(categorical_exposures, "globalcatmod")
  }
  cat("Sensitivity 'gba_binary': replaced globalbrainscore2 with globalcatmod (binary).\n")
}

# --- A.2 Metric-specific outlier exclusion (Methods 2.7.2) ---

if (length(subjects_to_exclude) > 0) {
  n_before <- nrow(df)
  df <- df[!df$ID %in% subjects_to_exclude, ]
  cat(sprintf("Excluded %d metric-specific outlier(s) by ID (%s). N = %d\n",
              n_before - nrow(df),
              paste(subjects_to_exclude, collapse = ", "), nrow(df)))
}

# --- A.3 Extract functional response ---

metric_cols <- grep(paste0("^", metric, "_"), colnames(df), value = TRUE)
cat(sprintf("\nFound %d density columns for %s\n", length(metric_cols), metric))
if (length(metric_cols) == 0) {
  stop(paste0("No columns matching '", metric, "_*'."))
}

densities_all <- as.numeric(gsub(paste0(metric, "_"), "", metric_cols))
keep <- densities_all >= density_min & densities_all <= density_max
densities <- densities_all[keep]
metric_cols <- metric_cols[keep]
n_density <- length(densities)

cat(sprintf("Density range: %.2f%% to %.2f%% (%d points)\n",
            min(densities), max(densities), n_density))

Y <- as.matrix(df[, metric_cols])
rownames(Y) <- df$ID
colnames(Y) <- densities

# --- A.4 Prepare predictors ---

all_vars <- unique(c("ID", forced_covariates, exposures))
predictor_df <- df %>% dplyr::select(all_of(all_vars))

# Apply GBA transform (only when continuous globalbrainscore2 is in the exposure
# set; for "gba_binary" branch the var has been swapped to globalcatmod and the
# runs driver sets gbs2_transform = "none")
if ("globalbrainscore2" %in% exposures && gbs2_transform != "none") {
  cat(sprintf("\nApplying '%s' transform to globalbrainscore2\n", gbs2_transform))
  cat(sprintf("  Before: range [%.1f, %.1f], median %.1f\n",
              min(predictor_df$globalbrainscore2, na.rm = TRUE),
              max(predictor_df$globalbrainscore2, na.rm = TRUE),
              median(predictor_df$globalbrainscore2, na.rm = TRUE)))
  if (gbs2_transform == "sqrt") {
    predictor_df$globalbrainscore2 <- sqrt(predictor_df$globalbrainscore2)
  } else if (gbs2_transform == "log1p") {
    predictor_df$globalbrainscore2 <- log1p(predictor_df$globalbrainscore2)
  } else if (gbs2_transform == "rank") {
    predictor_df$globalbrainscore2 <- rank(predictor_df$globalbrainscore2,
                                           na.last = "keep")
  }
  cat(sprintf("  After:  range [%.2f, %.2f], median %.2f\n",
              min(predictor_df$globalbrainscore2, na.rm = TRUE),
              max(predictor_df$globalbrainscore2, na.rm = TRUE),
              median(predictor_df$globalbrainscore2, na.rm = TRUE)))
}

# Complete cases on (ID + forced + exposures)
complete <- complete.cases(predictor_df)
Y_complete <- Y[complete, ]
predictor_df <- predictor_df[complete, ]
cat(sprintf("\nComplete cases: %d / %d\n", nrow(predictor_df), n_original))

# --- A.5 Standardize continuous variables ---

all_categorical <- unique(c(categorical_exposures, categorical_forced))
continuous_forced    <- setdiff(forced_covariates, all_categorical)
continuous_exposures <- setdiff(exposures, all_categorical)

standardization_params <- list()
predictor_df_scaled <- predictor_df

for (var in c(continuous_forced, continuous_exposures)) {
  if (var %in% colnames(predictor_df_scaled)) {
    m <- mean(predictor_df_scaled[[var]], na.rm = TRUE)
    s <- sd(predictor_df_scaled[[var]], na.rm = TRUE)
    standardization_params[[var]] <- list(mean = m, sd = s)
    predictor_df_scaled[[var]] <- (predictor_df_scaled[[var]] - m) / s
  }
}

cat(sprintf("\nStandardized:\n  Forced (continuous):    %s\n  Exposures (continuous): %s\n  Categorical:            %s\n",
            ifelse(length(continuous_forced) > 0,
                   paste(continuous_forced, collapse = ", "), "(none)"),
            ifelse(length(continuous_exposures) > 0,
                   paste(continuous_exposures, collapse = ", "), "(none)"),
            paste(intersect(all_categorical, c(forced_covariates, exposures)),
                  collapse = ", ")))

# Aliases used downstream (preserved from upstream pipeline naming)
Y_clean                   <- Y_complete
predictor_df_clean        <- predictor_df
predictor_df_scaled_clean <- predictor_df_scaled

n_final <- nrow(Y_clean)
cat(sprintf("\nFinal analysis sample: %d subjects\n", n_final))

# ==============================================================================
# SECTION B: STABILITY SELECTION VIA REPEATED SAMPLE SPLITTING
# ==============================================================================
#
# For each split b = 1, ..., n_splits:
#   1. Randomly split N subjects into selection (60%) / holdout (40%)
#   2. On selection set only:
#      a. PFFR-residualize Y on forced covariates
#      b. FWL-residualize each exposure on forced covariates
#      c. FPCA on Y-residuals; retain PCs by PVE threshold (with min/max)
#      d. Inverse-variance weights on PC scores
#      e. Cross-validated group LASSO on stacked design matrix; record
#         exposures with non-zero group L2 norm at lambda.min
#
# This eliminates selection-set / inference-set data leakage: FPCA basis
# and FWL residualization are fit only on the selection subset.
#
# ==============================================================================

cat("\n\n========================================\n")
cat("SECTION B: STABILITY SELECTION\n")
cat(sprintf("  %d splits, %.0f%%/%.0f%%, threshold %.0f%%\n",
            n_splits, selection_frac * 100, (1 - selection_frac) * 100,
            stability_threshold * 100))
cat("========================================\n\n")

n_exp <- length(exposures)
selection_matrix <- matrix(FALSE, nrow = n_splits, ncol = n_exp,
                           dimnames = list(NULL, exposures))

split_status <- character(n_splits)

# --- Helper: run one split ---

run_one_split <- function(split_idx, sel_idx, Y_clean,
                          predictor_df_scaled_clean, densities,
                          forced_covariates, exposures, pffr_k_basis,
                          fpca_pve_threshold, fpca_min_npc, fpca_max_npc,
                          n_cv_folds, pffr_family_final) {
  
  Y_sel    <- Y_clean[sel_idx, ]
  pred_sel <- predictor_df_scaled_clean[sel_idx, ]
  n_sel    <- nrow(Y_sel)
  
  # B.2: PFFR-residualize Y on forced covariates (selection set)
  pffr_data_sel <- pred_sel
  pffr_data_sel$Y <- Y_sel
  forced_formula <- as.formula(
    paste("Y ~", paste(forced_covariates, collapse = " + ")))
  
  pffr_sel <- tryCatch(
    pffr(forced_formula, yind = densities, data = pffr_data_sel,
         bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
         family = pffr_family_final),
    error = function(e) NULL
  )
  if (is.null(pffr_sel)) {
    return(list(retained = rep(FALSE, length(exposures)),
                status = "pffr_failed"))
  }
  
  Y_resid_sel <- residuals(pffr_sel)
  colnames(Y_resid_sel) <- densities
  
  # B.3: FWL-residualize each exposure on forced covariates (linear projection)
  X_exposures_raw <- as.matrix(pred_sel[, exposures, drop = FALSE])
  X_forced <- model.matrix(~ .,
                           data = pred_sel[, forced_covariates, drop = FALSE])[, -1, drop = FALSE]
  
  X_exposures_resid <- matrix(NA, nrow = n_sel, ncol = length(exposures))
  colnames(X_exposures_resid) <- exposures
  for (j in seq_along(exposures)) {
    fit_j <- lm(X_exposures_raw[, j] ~ X_forced)
    X_exposures_resid[, j] <- residuals(fit_j)
  }
  
  # B.4: FPCA on Y-residuals; retain PCs by PVE threshold with min/max bounds
  fpca_sel <- tryCatch(
    fpca.sc(Y = Y_resid_sel, argvals = densities,
            pve = fpca_pve_threshold, npc = fpca_max_npc, var = TRUE),
    error = function(e) NULL
  )
  if (is.null(fpca_sel)) {
    return(list(retained = rep(FALSE, length(exposures)),
                status = "fpca_failed"))
  }
  
  pve_individual <- fpca_sel$evalues / sum(fpca_sel$evalues) * 100
  pve_cumulative <- cumsum(pve_individual)
  n_pcs_pve <- which(pve_cumulative >= fpca_pve_threshold * 100)[1]
  if (is.na(n_pcs_pve)) n_pcs_pve <- fpca_sel$npc
  n_pcs <- max(fpca_min_npc, min(n_pcs_pve, fpca_max_npc, fpca_sel$npc))
  
  pc_scores_sel <- fpca_sel$scores[, 1:n_pcs, drop = FALSE]
  
  # B.5: inverse-variance PC weighting + group LASSO
  n_subj_sel  <- nrow(X_exposures_resid)
  n_exp_local <- ncol(X_exposures_resid)
  
  Y_stacked <- as.vector(pc_scores_sel)
  
  X_stacked <- matrix(0, nrow = n_subj_sel * n_pcs,
                      ncol = n_exp_local * n_pcs)
  for (k in 1:n_pcs) {
    row_idx <- ((k - 1) * n_subj_sel + 1):(k * n_subj_sel)
    col_idx <- ((1:n_exp_local) - 1) * n_pcs + k
    X_stacked[row_idx, col_idx] <- X_exposures_resid
  }
  
  group_labels <- rep(1:n_exp_local, each = n_pcs)
  
  # Inverse-variance PC weighting: PC noise variance = sum over densities of
  # (eigenfunction^2 * pointwise residual variance). Weights = 1/sqrt(noise_var),
  # normalized to mean 1 so the overall lambda scale is preserved.
  col_vars <- apply(Y_resid_sel, 2, stats::var)
  pc_noise_var <- vapply(1:n_pcs, function(k) {
    sum(fpca_sel$efunctions[, k]^2 * col_vars)
  }, numeric(1))
  pc_weights <- 1 / sqrt(pc_noise_var)
  pc_weights <- pc_weights / mean(pc_weights)
  
  for (k in 1:n_pcs) {
    row_idx <- ((k - 1) * n_subj_sel + 1):(k * n_subj_sel)
    Y_stacked[row_idx] <- Y_stacked[row_idx] * pc_weights[k]
    X_stacked[row_idx, ] <- X_stacked[row_idx, ] * pc_weights[k]
  }
  
  cv_fit <- tryCatch(
    cv.grpreg(X = X_stacked, y = Y_stacked, group = group_labels,
              penalty = "grLasso", nfolds = n_cv_folds, seed = split_idx),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) {
    return(list(retained = rep(FALSE, length(exposures)),
                status = "lasso_failed"))
  }
  
  coefs_sel <- coef(cv_fit, lambda = cv_fit$lambda.min)[-1]  # drop intercept
  retained <- vapply(1:n_exp_local, function(j) {
    any(abs(coefs_sel[group_labels == j]) > 1e-10)
  }, logical(1))
  
  list(retained = retained, status = "success")
}

# --- Run all splits ---

cat(sprintf("Running %d sample splits...\n", n_splits))
cat("Progress: ")

set.seed(master_seed)
split_seeds <- sample.int(1e7, n_splits)

for (b in 1:n_splits) {
  if (b %% 10 == 0) cat(sprintf("%d ", b))
  
  set.seed(split_seeds[b])
  n_sel   <- round(n_final * selection_frac)
  sel_idx <- sample(1:n_final, n_sel, replace = FALSE)
  
  result <- run_one_split(
    split_idx = b, sel_idx = sel_idx,
    Y_clean = Y_clean,
    predictor_df_scaled_clean = predictor_df_scaled_clean,
    densities = densities,
    forced_covariates = forced_covariates,
    exposures = exposures,
    pffr_k_basis = pffr_k_basis,
    fpca_pve_threshold = fpca_pve_threshold,
    fpca_min_npc = fpca_min_npc, fpca_max_npc = fpca_max_npc,
    n_cv_folds = n_cv_folds,
    pffr_family_final = pffr_family_final
  )
  
  selection_matrix[b, ] <- result$retained
  split_status[b]       <- result$status
}
cat("\nDone.\n\n")

# ==============================================================================
# SECTION C: SELECTION FREQUENCIES AND SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("SECTION C: SELECTION FREQUENCIES\n")
cat("========================================\n\n")

n_successful <- sum(split_status == "success")
cat(sprintf("Successful splits: %d / %d\n\n", n_successful, n_splits))
if (n_successful < n_splits) {
  cat("Failure breakdown:\n")
  print(table(split_status))
  cat("\n")
}

successful_mask <- split_status == "success"
selection_freq  <- colMeans(selection_matrix[successful_mask, , drop = FALSE])
names(selection_freq) <- exposures

# --- Save selection matrix (consumed by Supp Note 7) ---
sel_out <- selection_matrix[successful_mask, , drop = FALSE] * 1
stopifnot(identical(colnames(sel_out), exposures))

saveRDS(sel_out, file.path(output_dir, paste0(metric, "_selection_matrix.rds")))
write.csv(sel_out,
          file.path(output_dir, paste0(metric, "_selection_matrix.csv")),
          row.names = FALSE)
cat(sprintf("Saved selection matrix:\n  %s_selection_matrix.{rds,csv}\n\n", metric))

# --- Selection frequency summary ---
cat("--- Selection Frequencies ---\n\n")
for (j in seq_along(exposures)) {
  status <- ifelse(selection_freq[j] >= stability_threshold,
                   "** STABLE **", "  unstable")
  cat(sprintf("  %-25s  %.1f%%  (%d/%d splits)  %s\n",
              exposures[j], selection_freq[j] * 100,
              sum(selection_matrix[successful_mask, j]), n_successful, status))
}

stable_exposures   <- exposures[selection_freq >= stability_threshold]
unstable_exposures <- exposures[selection_freq <  stability_threshold]

cat(sprintf("\n  Stability threshold: %.0f%%\n", stability_threshold * 100))
cat(sprintf("  Stably selected: %s\n",
            ifelse(length(stable_exposures) > 0,
                   paste(stable_exposures, collapse = ", "), "(none)")))

# --- Selection frequency CSV (Table 3) ---
selection_summary <- data.frame(
  Exposure            = exposures,
  Selection_Frequency = selection_freq,
  Pct                 = round(selection_freq * 100, 1),
  Stable              = selection_freq >= stability_threshold,
  stringsAsFactors    = FALSE
)
write.csv(selection_summary,
          file.path(output_dir,
                    paste0(metric, "_stability_selection_summary.csv")),
          row.names = FALSE)
cat(sprintf("\nSaved: %s_stability_selection_summary.csv\n", metric))

# --- Selection frequency barplot (Figures 3A / 4A / 5) ---
freq_df <- data.frame(
  Exposure  = factor(exposures,
                     levels = exposures[order(selection_freq, decreasing = TRUE)]),
  Frequency = selection_freq * 100
)
p_freq <- ggplot(freq_df,
                 aes(x = Exposure, y = Frequency,
                     fill = Frequency >= stability_threshold * 100)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = stability_threshold * 100,
             linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60"),
                    labels = c("Below threshold", "Above threshold"),
                    name = "") +
  labs(x = "Exposure", y = "Selection Frequency (%)",
       title = paste(metric, "Stability Selection Frequencies"),
       subtitle = sprintf("%d splits (%.0f%%/%.0f%%), threshold = %.0f%%",
                          n_successful, selection_frac * 100,
                          (1 - selection_frac) * 100,
                          stability_threshold * 100)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir,
                 paste0(metric, "_stability_selection_frequencies.png")),
       p_freq, width = 8, height = 6, dpi = 300)
cat(sprintf("Saved: %s_stability_selection_frequencies.png\n", metric))

# ==============================================================================
# SECTION D: FULL-SAMPLE REFIT (CONDITIONAL ON SELECTION)
# ==============================================================================
#
# If any exposure was stably selected, fit pffr() on the full sample with
# forced covariates plus stably-selected exposures. Bootstrap CIs and per-
# predictor LR tests follow the same conventions as Methods 2.7.3.
#
# If no exposure was stably selected, Section D is skipped.
#
# ==============================================================================

cat("\n\n========================================\n")
cat("SECTION D: FULL-SAMPLE REFIT\n")
cat("  (conditional on stability selection)\n")
cat("========================================\n\n")

if (length(stable_exposures) == 0) {
  
  cat("No exposures met the stability threshold.\n")
  cat("Skipping full-sample refit per Methods 2.7.4 convention.\n")
  pffr_final          <- NULL
  bootstrap_coefs     <- NULL
  predictor_tests     <- NULL
  lr_test_overall     <- NULL
  ve_full             <- NA
  ve_forced           <- NA
  ve_intercept        <- NA
  functional_ve       <- NA
  functional_ve_forced <- NA
  functional_ve_increment <- NA
  overall_p           <- NA
  overall_test_stat   <- NA
  overall_significant <- NA
  lr_p_vs_forced      <- NA
  lr_stat_vs_forced   <- NA
  final_predictors    <- forced_covariates
  
} else {
  
  final_predictors <- c(forced_covariates, stable_exposures)
  cat(sprintf("Final model predictors (%d forced + %d selected):\n  %s\n\n",
              length(forced_covariates), length(stable_exposures),
              paste(final_predictors, collapse = ", ")))
  
  # ---- D.1 Fit final model, forced-only reference, intercept-only reference ----
  pffr_data_final <- predictor_df_scaled_clean[, c("ID", final_predictors)]
  pffr_data_final$Y <- Y_clean
  
  final_formula  <- as.formula(paste("Y ~",
                                     paste(final_predictors, collapse = " + ")))
  forced_formula <- as.formula(paste("Y ~",
                                     paste(forced_covariates, collapse = " + ")))
  
  cat(sprintf("Fitting pffr() final model (family = %s)...\n",
              if (inherits(pffr_family_final, "general.family"))
                "scat" else pffr_family_final$family))
  
  pffr_final <- pffr(
    final_formula, yind = densities, data = pffr_data_final,
    bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
    family = pffr_family_final
  )
  
  pffr_forced_full <- pffr(
    forced_formula, yind = densities, data = pffr_data_final,
    bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
    family = pffr_family_final
  )
  
  pffr_intercept <- pffr(
    Y ~ 1, yind = densities, data = pffr_data_final,
    bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
    family = pffr_family_final
  )
  
  cat("Final, forced-only, and intercept-only models fit.\n\n")
  print(summary(pffr_final))
  
  # ---- D.2 Variance explained (deviance explained for non-gaussian) ----
  get_var_explained <- function(model) {
    s <- summary(model)
    if (lr_test_type == "F") s$r.sq else s$dev.expl
  }
  
  ve_full     <- get_var_explained(pffr_final)
  ve_forced   <- get_var_explained(pffr_forced_full)
  ve_intercept <- get_var_explained(pffr_intercept)
  
  functional_ve           <- (ve_full   - ve_intercept) / (1 - ve_intercept)
  functional_ve_forced    <- (ve_forced - ve_intercept) / (1 - ve_intercept)
  functional_ve_increment <- functional_ve - functional_ve_forced
  ve_label <- if (lr_test_type == "F") "R^2" else "Deviance explained"
  
  cat(sprintf("\n%s (full):                   %.4f (%.2f%%)\n",
              ve_label, ve_full,    ve_full    * 100))
  cat(sprintf("%s (forced):                 %.4f (%.2f%%)\n",
              ve_label, ve_forced,  ve_forced  * 100))
  cat(sprintf("Functional %s (full):         %.4f (%.2f%%)\n",
              ve_label, functional_ve, functional_ve * 100))
  cat(sprintf("Functional %s (forced):       %.4f (%.2f%%)\n",
              ve_label, functional_ve_forced, functional_ve_forced * 100))
  cat(sprintf("Functional %s increment:      %.4f (%.2f%%)\n\n",
              ve_label, functional_ve_increment, functional_ve_increment * 100))
  
  # ---- D.3 Overall LR test (varying vs. constant coefficients) ----
  constant_parts   <- paste0("c(", final_predictors, ")", collapse = " + ")
  constant_formula <- as.formula(paste("Y ~ 1 +", constant_parts))
  
  pffr_constant <- pffr(
    constant_formula, yind = densities, data = pffr_data_final,
    bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
    family = pffr_family_final
  )
  
  mtest_final    <- pffr_final
  mtest_constant <- pffr_constant
  class(mtest_final)    <- class(mtest_final)[-1]
  class(mtest_constant) <- class(mtest_constant)[-1]
  
  lr_test_overall <- anova(mtest_constant, mtest_final, test = lr_test_type)
  
  cat("--- Overall LR Test (varying vs constant coefficients) ---\n")
  if (lr_test_type == "F") {
    overall_test_stat <- lr_test_overall$F[2]
    overall_p         <- lr_test_overall$`Pr(>F)`[2]
    cat(sprintf("  F = %.2f, p = %.2e\n", overall_test_stat, overall_p))
  } else {
    overall_test_stat <- lr_test_overall$Deviance[2]
    overall_p         <- lr_test_overall$`Pr(>Chi)`[2]
    cat(sprintf("  Chi^2 = %.4f, p = %.2e\n", overall_test_stat, overall_p))
  }
  overall_significant <- overall_p < p_threshold
  cat(sprintf("  Significant at p < %.3f: %s\n\n", p_threshold,
              ifelse(overall_significant, "YES", "NO")))
  
  # ---- D.4 LR test: full vs forced-only ----
  mtest_forced <- pffr_forced_full
  class(mtest_forced) <- class(mtest_forced)[-1]
  
  lr_test_vs_forced <- anova(mtest_forced, mtest_final, test = lr_test_type)
  cat("--- LR Test (selected exposures vs forced-only) ---\n")
  if (lr_test_type == "F") {
    lr_stat_vs_forced <- lr_test_vs_forced$F[2]
    lr_p_vs_forced    <- lr_test_vs_forced$`Pr(>F)`[2]
    cat(sprintf("  F = %.2f, p = %.2e\n\n", lr_stat_vs_forced, lr_p_vs_forced))
  } else {
    lr_stat_vs_forced <- lr_test_vs_forced$Deviance[2]
    lr_p_vs_forced    <- lr_test_vs_forced$`Pr(>Chi)`[2]
    cat(sprintf("  Chi^2 = %.2f, p = %.2e\n\n",
                lr_stat_vs_forced, lr_p_vs_forced))
  }
  
  # ---- D.5 Per-predictor LR tests ----
  cat("--- Per-Predictor LR Tests ---\n\n")
  predictor_tests <- data.frame(
    Predictor       = character(),
    Test_Statistic  = numeric(),
    p_value         = numeric(),
    functional_ve   = numeric(),
    semi_partial_r  = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pred in final_predictors) {
    other_preds <- setdiff(final_predictors, pred)
    reduced_formula <- as.formula(
      paste("Y ~", paste(other_preds, collapse = " + ")))
    
    pffr_reduced <- tryCatch(
      pffr(reduced_formula, yind = densities, data = pffr_data_final,
           bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
           family = pffr_family_final),
      error = function(e) NULL
    )
    if (is.null(pffr_reduced)) {
      cat(sprintf("  %s: Reduced model fit failed; skipping.\n", pred)); next
    }
    mtest_reduced <- pffr_reduced
    class(mtest_reduced) <- class(mtest_reduced)[-1]
    
    lr_test <- tryCatch(
      anova(mtest_reduced, mtest_final, test = lr_test_type),
      error = function(e) NULL
    )
    if (is.null(lr_test)) {
      cat(sprintf("  %s: LR test failed; skipping.\n", pred)); next
    }
    
    single_formula <- as.formula(paste("Y ~ 1 +", pred))
    pffr_single <- tryCatch(
      pffr(single_formula, yind = densities, data = pffr_data_final,
           bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
           family = pffr_family_final),
      error = function(e) NULL
    )
    func_ve_pred <- if (!is.null(pffr_single))
      (get_var_explained(pffr_single) - ve_intercept) / (1 - ve_intercept)
    else NA
    
    ve_reduced     <- get_var_explained(pffr_reduced)
    semi_partial_r <- sqrt(abs(ve_full - ve_reduced))
    
    pred_stat <- if (lr_test_type == "F") lr_test$F[2]
    else                     lr_test$Deviance[2]
    pred_p    <- if (lr_test_type == "F") lr_test$`Pr(>F)`[2]
    else                     lr_test$`Pr(>Chi)`[2]
    
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
                pred, if (lr_test_type == "F") "F" else "Chi^2",
                pred_stat, pred_p, ve_label, func_ve_pred,
                semi_partial_r, sig_marker))
  }
  cat(sprintf("\n  Significance: *** p < %.3f, ** p < 0.01, * p < 0.05\n",
              p_threshold))
  
  write.csv(predictor_tests,
            file.path(output_dir, paste0(metric, "_predictor_tests.csv")),
            row.names = FALSE)
  cat(sprintf("\n  Saved: %s_predictor_tests.csv\n", metric))
  
  # ---- D.6 Bootstrap CIs ----
  cat(sprintf("\nRunning bootstrap (%d iterations)...\n", n_bootstrap))
  n_cores <- max(1, detectCores() - 1)
  cat(sprintf("Using %d cores for parallel processing.\n", n_cores))
  
  bootstrap_coefs <- tryCatch({
    coefboot.pffr(pffr_final, B = n_bootstrap,
                  ncpus = n_cores, parallel = "multicore")
  }, error = function(e) {
    cat(sprintf("coefboot.pffr error: %s\n", e$message)); NULL
  })
  if (!is.null(bootstrap_coefs)) cat("Bootstrap complete.\n")
  
  # ---- D.7 Bootstrap CI plots ----
  if (!is.null(bootstrap_coefs)) {
    
    getPlotObject <- function(model) {
      ff <- tempfile(); svg(filename = ff)
      plotObject <- plot(model)
      dev.off(); unlink(ff)
      plotObject
    }
    
    getCIsList <- function(coefboot_bs) {
      smList <- coefboot_bs$smterms
      lapply(smList, function(s) {
        list(x = s[[2]], y = s[[1]],
             ci_lower = s$`2.5%`, ci_upper = s$`97.5%`)
      })
    }
    
    CIs_list  <- getCIsList(bootstrap_coefs)
    variable_names_list <- c("Intercept", final_predictors)
    plotObject <- getPlotObject(pffr_final)
    
    for (i in seq_along(plotObject)) {
      pv <- plotObject[[i]]
      plot_df <- data.frame(
        x = pv$x, y = pv$fit,
        ci_lower = CIs_list[[i]]$ci_lower,
        ci_upper = CIs_list[[i]]$ci_upper
      )
      gg <- ggplot(plot_df, aes(x = x, y = y)) +
        geom_hline(yintercept = 0, lty = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    alpha = 0.3, fill = "steelblue") +
        geom_line(linewidth = 1, color = "darkblue") +
        ylab(expression(beta(d))) + xlab("Density (%)") +
        ggtitle(paste0(variable_names_list[[i]],
                       " Effect on ", metric, " (95% Bootstrap CIs)")) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file.path(output_dir,
                       paste0(metric, "_final_beta_",
                              gsub("[^a-zA-Z0-9]", "_", variable_names_list[[i]]),
                              ".png")),
             gg, width = 8, height = 5, dpi = 300)
    }
    cat("Saved bootstrap CI plots.\n")
  }
  
  # ---- D.8 Model diagnostics (Methods 2.7.3 / 2.7.4: k-basis, residuals) ----
  png(file.path(output_dir, paste0(metric, "_pffr_diagnostics.png")),
      width = 12, height = 10, units = "in", res = 300)
  pffr.check(pffr_final)
  dev.off()
  cat(sprintf("Saved: %s_pffr_diagnostics.png\n", metric))
  
  Y_fitted        <- fitted(pffr_final)
  residuals_pffr  <- Y_clean - Y_fitted
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
  ggsave(file.path(output_dir,
                   paste0(metric, "_final_pffr_residuals.png")),
         p_resid, width = 10, height = 6, dpi = 300)
  cat(sprintf("Saved: %s_final_pffr_residuals.png\n", metric))
}

# ==============================================================================
# SECTION E: SAVE RESULTS
# ==============================================================================

cat("\n\n========================================\n")
cat("SECTION E: SAVE RESULTS\n")
cat("========================================\n\n")

ve_label <- if (lr_test_type == "F") "R^2" else "Deviance explained"

model_summary <- data.frame(
  Metric                = metric,
  Sensitivity_Branch    = sensitivity_branch,
  Density_Range         = paste0(density_min, "-", density_max, "%"),
  N_Final               = n_final,
  Forced_Covariates     = paste(forced_covariates, collapse = ", "),
  N_Exposures_Tested    = n_exp,
  N_Splits              = n_splits,
  N_Successful_Splits   = n_successful,
  Stability_Threshold   = stability_threshold,
  N_Stable_Exposures    = length(stable_exposures),
  Stable_Exposures      = paste(stable_exposures,   collapse = ", "),
  Unstable_Exposures    = paste(unstable_exposures, collapse = ", "),
  Family                = if (inherits(pffr_family_final, "general.family"))
    "scat" else pffr_family_final$family,
  LR_Test_Type          = lr_test_type,
  VE_Label              = ve_label,
  VE_Full               = ve_full,
  VE_Forced_Only        = ve_forced,
  Functional_VE         = functional_ve,
  Functional_VE_Forced  = functional_ve_forced,
  Functional_VE_Increment = functional_ve_increment,
  Overall_Test_Stat     = overall_test_stat,
  Overall_p             = overall_p,
  Overall_Significant   = overall_significant,
  LR_vs_Forced_Stat     = lr_stat_vs_forced,
  LR_vs_Forced_p        = lr_p_vs_forced,
  GBA_Transform         = gbs2_transform,
  Bootstrap_N           = n_bootstrap,
  Bootstrap_Completed   = !is.null(bootstrap_coefs),
  stringsAsFactors      = FALSE
)
write.csv(model_summary,
          file.path(output_dir, paste0(metric, "_model_summary.csv")),
          row.names = FALSE)
cat(sprintf("Saved: %s_model_summary.csv\n", metric))

save(
  metric, sensitivity_branch, densities, density_min, density_max,
  forced_covariates, exposures, stable_exposures, unstable_exposures,
  pffr_family_final, lr_test_type, ve_label,
  Y_clean, predictor_df_clean, predictor_df_scaled_clean,
  standardization_params,
  selection_matrix, selection_freq, split_status,
  selection_summary, model_summary,
  pffr_final, bootstrap_coefs, predictor_tests, lr_test_overall,
  ve_full, ve_forced, ve_intercept,
  functional_ve, functional_ve_forced, functional_ve_increment,
  overall_test_stat, overall_p, overall_significant,
  lr_stat_vs_forced, lr_p_vs_forced,
  file = file.path(output_dir, paste0(metric, "_stabsel_results.RData"))
)
cat(sprintf("Saved: %s_stabsel_results.RData\n", metric))

cat("\n\nSession Information:\n")
print(sessionInfo())

# ==============================================================================
# RUN SUMMARY
# ==============================================================================

cat("\n\n========================================\n")
cat("RUN SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Metric:              %s\n", metric))
cat(sprintf("Sensitivity branch:  %s\n", sensitivity_branch))
cat(sprintf("Sample size:         %d\n", n_final))
cat(sprintf("Stably selected:     %s\n",
            ifelse(length(stable_exposures) > 0,
                   paste(stable_exposures, collapse = ", "), "(none)")))

if (length(stable_exposures) > 0) {
  cat(sprintf("\n%s (full):       %.4f\n", ve_label, ve_full))
  cat(sprintf("Functional %s:    %.4f (%.2f%%)\n",
              ve_label, functional_ve, functional_ve * 100))
  cat(sprintf("Functional %s increment from selected exposures: %.4f (%.2f%%)\n",
              ve_label, functional_ve_increment, functional_ve_increment * 100))
  cat(sprintf("\nOverall LR test: %s = %.2f, p = %.2e (significant at p < %.3f: %s)\n",
              if (lr_test_type == "F") "F" else "Chi^2",
              overall_test_stat, overall_p, p_threshold,
              ifelse(overall_significant, "YES", "NO")))
  cat(sprintf("LR test vs forced-only: %s = %.2f, p = %.2e\n",
              if (lr_test_type == "F") "F" else "Chi^2",
              lr_stat_vs_forced, lr_p_vs_forced))
}

cat(sprintf("\nOutputs in:\n%s\n", output_dir))