################################################################################
#
#   note7_stability_selection_sensitivity.R  -  Supplementary Note 7:
#   robustness of stability selection to its parameter settings
#
#   Re-runs ONLY the stability-selection step (Sections A-C of
#   code/03_statistical_analysis/FOUR_fda_stability_selection_pipeline.R)
#   under alternative settings.
#   No full-sample refits, no bootstrap.
#
#   The data preparation and the per-split procedure are copied unchanged from
#   the pipeline. The only addition is that each split also records the number
#   of FPCs retained, so the PVE sensitivity can be described.
#
#   CONFIGURATIONS (one setting varied at a time; all others at manuscript values)
#     Block 1  more iterations : 200 splits, 60% subsample, PVE 0.995
#     Block 2  subsample       : 200 splits, subsample 0.50 / 0.70 / 0.80, PVE 0.995
#     Block 3  FPC retention   : 200 splits, 60% subsample, PVE 0.90 / 0.95 / 0.99
#     Each block x 4 metrics (strength, normalized GE, normalized ACC, SW).
#     Selection-frequency thresholds are swept afterwards (--summarize); they
#     need no new runs.
#
#   USAGE (Terminal, from the repository root)
#     Rscript code/05_supplemental/note7_stability_selection_sensitivity.R --list        # configuration table
#     Rscript code/05_supplemental/note7_stability_selection_sensitivity.R all 8         # every config, 8 cores
#     Rscript code/05_supplemental/note7_stability_selection_sensitivity.R 1-4 8         # configs 1 to 4
#     Rscript code/05_supplemental/note7_stability_selection_sensitivity.R 6 8           # config 6 only
#     Rscript code/05_supplemental/note7_stability_selection_sensitivity.R --summarize   # summary tables
#   Finished configs (summary.csv present) are skipped, so an interrupted run
#   can be restarted with the same command.
#   Results are identical for any number of cores (each split is seeded).
#
#   OUTPUT  (OUT_ROOT = results/stability_selection_sensitivity)
#     {OUT_ROOT}/{metric_folder}/{config}/
#        selection_matrix.csv   splits x exposures (1 = retained), successful splits
#        split_log.csv          split, status, n_pcs
#        summary.csv            exposure, selection frequency
#     {OUT_ROOT}/sensitivity_selection_frequencies.csv   all configs, long format
#     {OUT_ROOT}/sensitivity_threshold_sweep.csv         stable exposures by threshold
#
################################################################################

# ------------------------------------------------------------------ CONFIG
# Repo root is located automatically so a fresh clone runs without edits.
# Order: (0) honor a `repo_root` already set in the global env; (1) derive this
# script's own location and walk upward; (2) walk upward from the working dir.
# A folder is the repo root if it contains both `data/analysis_ready/` and `code/`.
# Manual override if auto-detection ever fails:
#   repo_root <- "/full/path/to/repository"   # then source() this script
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
  repo_root <- NULL
  this_script <- tryCatch({
    sf <- sys.frames(); paths <- character(0)
    for (fr in sf) {
      ofile <- tryCatch(get("ofile", envir = fr, inherits = FALSE), error = function(e) NULL)
      if (!is.null(ofile) && is.character(ofile) && nzchar(ofile)) paths <- c(paths, ofile)
    }
    if (length(paths) > 0) paths[length(paths)] else NULL
  }, error = function(e) NULL)
  if (!is.null(this_script)) repo_root <- .find_repo_root_from_path(dirname(this_script))
  if (is.null(repo_root))    repo_root <- .find_repo_root_from_path(getwd())
  if (is.null(repo_root))
    stop("Could not locate repo root. Set it manually before source()-ing:\n",
         "  repo_root <- \"/full/path/to/repository\"")
}
cat(sprintf("Repo root: %s\n", repo_root))

DATA_PATH <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
OUT_ROOT  <- file.path(repo_root, "results/stability_selection_sensitivity")
MAIN_ROOT <- file.path(repo_root, "results/stability_selection")   # existing 100-split runs (read by --summarize)

FORCED      <- c("eTIV", "sex", "sriskscore", "age_at_5y_mri", "Rel_Motion")
EXPOSURES   <- c("bpd2", "bw_z_new", "ga", "globalbrainscore2", "anyrop", "sepsis2", "anyivh", "hydrocephalus_dc")
CAT_EXP     <- c("bpd2", "anyrop", "sepsis2", "anyivh", "hydrocephalus_dc")
CAT_FORCED  <- c("sex")

METRICS <- list(
  list(metric = "str",               folder = "strength", family = "gaussian", exclude = c()),
  list(metric = "rand_norm_wei_GE",  folder = "GE_norm",  family = "gaussian", exclude = c(713, 659)),
  list(metric = "rand_norm_wei_ACC", folder = "ACC_norm", family = "gaussian", exclude = c()),
  list(metric = "rand_norm_wei_SW",  folder = "SW",       family = "scat",     exclude = c())
)

SETTINGS <- list(
  list(name = "splits200",           n_splits = 200, sel_frac = 0.60, pve = 0.995),
  list(name = "splits200_sub0.50",   n_splits = 200, sel_frac = 0.50, pve = 0.995),
  list(name = "splits200_sub0.70",   n_splits = 200, sel_frac = 0.70, pve = 0.995),
  list(name = "splits200_sub0.80",   n_splits = 200, sel_frac = 0.80, pve = 0.995),
  list(name = "splits200_pve0.90",   n_splits = 200, sel_frac = 0.60, pve = 0.90),
  list(name = "splits200_pve0.95",   n_splits = 200, sel_frac = 0.60, pve = 0.95),
  list(name = "splits200_pve0.99",   n_splits = 200, sel_frac = 0.60, pve = 0.99)
)

THRESHOLDS <- c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)

# Fixed at manuscript values (as in the pipeline)
DENSITY_MIN <- 11; DENSITY_MAX <- 100
FPCA_MIN_NPC <- 3; FPCA_MAX_NPC <- 10
PFFR_K <- 20; N_CV_FOLDS <- 10; MASTER_SEED <- 42

# ------------------------------------------------------------------ CONFIG TABLE
configs <- list()
for (s in SETTINGS) for (m in METRICS) configs[[length(configs) + 1]] <- c(m, s)
cfg_table <- do.call(rbind, lapply(seq_along(configs), function(i) {
  x <- configs[[i]]
  data.frame(index = i, metric = x$folder, config = x$name, n_splits = x$n_splits,
             subsample = x$sel_frac, pve = x$pve, stringsAsFactors = FALSE)
}))
out_dir_for <- function(x) file.path(OUT_ROOT, x$folder, x$name)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript note7_stability_selection_sensitivity.R --list | --summarize | all|<i>|<i-j> [cores]")

if (args[1] == "--list") {
  cfg_table$done <- vapply(configs, function(x) file.exists(file.path(out_dir_for(x), "summary.csv")), logical(1))
  print(cfg_table, row.names = FALSE)
  quit(status = 0)
}

# ------------------------------------------------------------------ SUMMARIZE
if (args[1] == "--summarize") {
  rows <- list()
  add_matrix <- function(path, folder, name, n_splits, sel_frac, pve) {
    if (!file.exists(path)) return(invisible(NULL))
    sm <- read.csv(path, check.names = FALSE)
    rows[[length(rows) + 1]] <<- data.frame(
      metric = folder, config = name, n_splits = n_splits, n_successful = nrow(sm),
      subsample = sel_frac, pve = pve, exposure = colnames(sm),
      selection_frequency = round(colMeans(sm), 3), stringsAsFactors = FALSE)
  }
  metric_prefix <- setNames(vapply(METRICS, `[[`, "", "metric"), vapply(METRICS, `[[`, "", "folder"))
  for (m in METRICS)   # existing 100-split runs at manuscript settings
    add_matrix(file.path(MAIN_ROOT, m$folder, paste0(m$metric, "_selection_matrix.csv")),
               m$folder, "splits100 (main)", 100, 0.60, 0.995)
  for (x in configs)
    add_matrix(file.path(out_dir_for(x), "selection_matrix.csv"),
               x$folder, x$name, x$n_splits, x$sel_frac, x$pve)
  if (length(rows) == 0) stop("No results found yet.")
  freq <- do.call(rbind, rows); rownames(freq) <- NULL
  write.csv(freq, file.path(OUT_ROOT, "sensitivity_selection_frequencies.csv"), row.names = FALSE)

  sweep <- do.call(rbind, lapply(split(freq, list(freq$metric, freq$config), drop = TRUE), function(d) {
    do.call(rbind, lapply(THRESHOLDS, function(t) data.frame(
      metric = d$metric[1], config = d$config[1], threshold = t,
      stable_exposures = { s <- d$exposure[d$selection_frequency >= t]; if (length(s)) paste(s, collapse = ", ") else "(none)" },
      stringsAsFactors = FALSE)))
  }))
  rownames(sweep) <- NULL
  sweep <- sweep[order(sweep$metric, sweep$config, sweep$threshold), ]
  write.csv(sweep, file.path(OUT_ROOT, "sensitivity_threshold_sweep.csv"), row.names = FALSE)
  cat(sprintf("Wrote %s\n      %s\n", file.path(OUT_ROOT, "sensitivity_selection_frequencies.csv"),
              file.path(OUT_ROOT, "sensitivity_threshold_sweep.csv")))
  wide <- reshape(freq[, c("metric", "config", "exposure", "selection_frequency")],
                  idvar = c("metric", "exposure"), timevar = "config", direction = "wide")
  names(wide) <- sub("^selection_frequency\\.", "", names(wide))
  print(wide[order(wide$metric, wide$exposure), ], row.names = FALSE)
  quit(status = 0)
}

# ------------------------------------------------------------------ WHICH CONFIGS
sel <- if (args[1] == "all") seq_along(configs) else if (grepl("^[0-9]+-[0-9]+$", args[1])) {
  ab <- as.integer(strsplit(args[1], "-")[[1]]); ab[1]:ab[2]
} else as.integer(args[1])
if (any(is.na(sel)) || any(sel < 1 | sel > length(configs))) stop("Config index out of range 1-", length(configs))
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

suppressPackageStartupMessages({
  library(refund); library(mgcv); library(readxl); library(dplyr); library(grpreg); library(parallel)
})

# ------------------------------------------------------------------ DATA PREP (pipeline Section A, main branch)
prep_data <- function(metric, subjects_to_exclude) {
  df <- read_excel(DATA_PATH)
  if (length(subjects_to_exclude) > 0) df <- df[!df$ID %in% subjects_to_exclude, ]
  metric_cols <- grep(paste0("^", metric, "_"), colnames(df), value = TRUE)
  densities_all <- as.numeric(gsub(paste0(metric, "_"), "", metric_cols))
  keep <- densities_all >= DENSITY_MIN & densities_all <= DENSITY_MAX
  densities <- densities_all[keep]; metric_cols <- metric_cols[keep]
  Y <- as.matrix(df[, metric_cols]); rownames(Y) <- df$ID; colnames(Y) <- densities
  predictor_df <- df %>% dplyr::select(all_of(unique(c("ID", FORCED, EXPOSURES))))
  predictor_df$globalbrainscore2 <- sqrt(predictor_df$globalbrainscore2)        # gbs2_transform = "sqrt"
  complete <- complete.cases(predictor_df)
  Y <- Y[complete, ]; predictor_df <- predictor_df[complete, ]
  all_cat <- unique(c(CAT_EXP, CAT_FORCED))
  for (v in c(setdiff(FORCED, all_cat), setdiff(EXPOSURES, all_cat))) {
    predictor_df[[v]] <- (predictor_df[[v]] - mean(predictor_df[[v]], na.rm = TRUE)) / sd(predictor_df[[v]], na.rm = TRUE)
  }
  list(Y = Y, X = predictor_df, densities = densities)
}

# ------------------------------------------------------------------ ONE SPLIT (pipeline run_one_split, + n_pcs)
run_one_split <- function(split_idx, sel_idx, Y_clean, predictor_df_scaled_clean, densities,
                          forced_covariates, exposures, pffr_k_basis, fpca_pve_threshold,
                          fpca_min_npc, fpca_max_npc, n_cv_folds, pffr_family_final) {
  fail <- function(st) list(retained = rep(FALSE, length(exposures)), status = st, n_pcs = NA_integer_)
  Y_sel <- Y_clean[sel_idx, ]; pred_sel <- predictor_df_scaled_clean[sel_idx, ]; n_sel <- nrow(Y_sel)

  pffr_data_sel <- pred_sel; pffr_data_sel$Y <- Y_sel
  forced_formula <- as.formula(paste("Y ~", paste(forced_covariates, collapse = " + ")))
  pffr_sel <- tryCatch(pffr(forced_formula, yind = densities, data = pffr_data_sel,
                            bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
                            family = pffr_family_final), error = function(e) NULL)
  if (is.null(pffr_sel)) return(fail("pffr_failed"))
  Y_resid_sel <- residuals(pffr_sel); colnames(Y_resid_sel) <- densities

  X_exposures_raw <- as.matrix(pred_sel[, exposures, drop = FALSE])
  X_forced <- model.matrix(~ ., data = pred_sel[, forced_covariates, drop = FALSE])[, -1, drop = FALSE]
  X_exposures_resid <- matrix(NA, nrow = n_sel, ncol = length(exposures)); colnames(X_exposures_resid) <- exposures
  for (j in seq_along(exposures)) X_exposures_resid[, j] <- residuals(lm(X_exposures_raw[, j] ~ X_forced))

  fpca_sel <- tryCatch(fpca.sc(Y = Y_resid_sel, argvals = densities, pve = fpca_pve_threshold,
                               npc = fpca_max_npc, var = TRUE), error = function(e) NULL)
  if (is.null(fpca_sel)) return(fail("fpca_failed"))
  pve_cumulative <- cumsum(fpca_sel$evalues / sum(fpca_sel$evalues) * 100)
  n_pcs_pve <- which(pve_cumulative >= fpca_pve_threshold * 100)[1]
  if (is.na(n_pcs_pve)) n_pcs_pve <- fpca_sel$npc
  n_pcs <- max(fpca_min_npc, min(n_pcs_pve, fpca_max_npc, fpca_sel$npc))
  pc_scores_sel <- fpca_sel$scores[, 1:n_pcs, drop = FALSE]

  n_subj_sel <- nrow(X_exposures_resid); n_exp_local <- ncol(X_exposures_resid)
  Y_stacked <- as.vector(pc_scores_sel)
  X_stacked <- matrix(0, nrow = n_subj_sel * n_pcs, ncol = n_exp_local * n_pcs)
  for (k in 1:n_pcs) {
    row_idx <- ((k - 1) * n_subj_sel + 1):(k * n_subj_sel)
    X_stacked[row_idx, ((1:n_exp_local) - 1) * n_pcs + k] <- X_exposures_resid
  }
  group_labels <- rep(1:n_exp_local, each = n_pcs)
  col_vars <- apply(Y_resid_sel, 2, stats::var)
  pc_noise_var <- vapply(1:n_pcs, function(k) sum(fpca_sel$efunctions[, k]^2 * col_vars), numeric(1))
  pc_weights <- 1 / sqrt(pc_noise_var); pc_weights <- pc_weights / mean(pc_weights)
  for (k in 1:n_pcs) {
    row_idx <- ((k - 1) * n_subj_sel + 1):(k * n_subj_sel)
    Y_stacked[row_idx] <- Y_stacked[row_idx] * pc_weights[k]
    X_stacked[row_idx, ] <- X_stacked[row_idx, ] * pc_weights[k]
  }
  cv_fit <- tryCatch(cv.grpreg(X = X_stacked, y = Y_stacked, group = group_labels, penalty = "grLasso",
                               nfolds = n_cv_folds, seed = split_idx), error = function(e) NULL)
  if (is.null(cv_fit)) return(fail("lasso_failed"))
  coefs_sel <- coef(cv_fit, lambda = cv_fit$lambda.min)[-1]
  retained <- vapply(1:n_exp_local, function(j) any(abs(coefs_sel[group_labels == j]) > 1e-10), logical(1))
  list(retained = retained, status = "success", n_pcs = as.integer(n_pcs))
}

# ------------------------------------------------------------------ RUN
for (i in sel) {
  x <- configs[[i]]; od <- out_dir_for(x)
  if (file.exists(file.path(od, "summary.csv"))) { cat(sprintf("[%d] %s / %s already done, skipping\n", i, x$folder, x$name)); next }
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n[%d/%d] %s | %s | %d splits, subsample %.2f, PVE %.3f | %d cores | started %s\n",
              i, length(configs), x$folder, x$name, x$n_splits, x$sel_frac, x$pve, n_cores, format(Sys.time())))

  d <- prep_data(x$metric, x$exclude)
  fam <- if (x$family == "scat") scat() else gaussian()
  n_final <- nrow(d$Y); n_sel <- round(n_final * x$sel_frac)
  set.seed(MASTER_SEED); split_seeds <- sample.int(1e7, x$n_splits)
  cat(sprintf("  N = %d, selection subset = %d\n", n_final, n_sel))

  one <- function(b) {
    set.seed(split_seeds[b])
    sel_idx <- sample(1:n_final, n_sel, replace = FALSE)
    run_one_split(b, sel_idx, d$Y, d$X, d$densities, FORCED, EXPOSURES, PFFR_K, x$pve,
                  FPCA_MIN_NPC, FPCA_MAX_NPC, N_CV_FOLDS, fam)
  }
  res <- if (n_cores > 1) mclapply(seq_len(x$n_splits), one, mc.cores = n_cores, mc.preschedule = FALSE)
         else lapply(seq_len(x$n_splits), one)

  bad <- vapply(res, function(r) !is.list(r) || is.null(r$status), logical(1))
  if (any(bad)) res[bad] <- lapply(seq_len(sum(bad)), function(k) list(retained = rep(FALSE, length(EXPOSURES)), status = "worker_error", n_pcs = NA_integer_))
  status <- vapply(res, `[[`, "", "status")
  selm <- do.call(rbind, lapply(res, `[[`, "retained")) * 1; colnames(selm) <- EXPOSURES
  ok <- status == "success"
  write.csv(selm[ok, , drop = FALSE], file.path(od, "selection_matrix.csv"), row.names = FALSE)
  write.csv(data.frame(split = seq_len(x$n_splits), status = status,
                       n_pcs = vapply(res, function(r) as.integer(r$n_pcs), integer(1))),
            file.path(od, "split_log.csv"), row.names = FALSE)
  fr <- colMeans(selm[ok, , drop = FALSE])
  write.csv(data.frame(exposure = EXPOSURES, selection_frequency = round(fr, 3), n_successful = sum(ok)),
            file.path(od, "summary.csv"), row.names = FALSE)
  cat(sprintf("  successful splits %d/%d | finished %s\n", sum(ok), x$n_splits, format(Sys.time())))
  for (j in seq_along(EXPOSURES)) cat(sprintf("    %-18s %5.1f%%\n", EXPOSURES[j], 100 * fr[j]))
}
cat("\nDone. Run with --summarize to build the summary tables.\n")
