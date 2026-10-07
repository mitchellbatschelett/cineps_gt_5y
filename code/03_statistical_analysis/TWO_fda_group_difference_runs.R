################################################################################
#
#   TWO_fda_group_difference_runs.R  -  VPT vs FT group-difference models
#
#   36 runs = 6 metrics x 6 covariate configurations, each fit by fda_pipeline.R
#     configurations: 1 unadjusted, 2 +age, 3 +eTIV, 4 +sex, 5 +motion, 6 fully adjusted
#     family: scat() for SW, gaussian() for every other metric
#     motion covariate: Rel_Motion (eddy relative RMS)
#     z-scored: age_at_5y_mri, eTIV, Rel_Motion; not scaled: Group, sex
#     metric-specific exclusions: rand_norm_wei_GE c(713, 659); none for other metrics
#
#   INDEX TABLE (SW last)
#     1-6 strength | 7-12 GE_norm | 13-18 ACC_norm | 19-24 GE_raw | 25-30 ACC_raw | 31-36 SW (scat)
#     within each block of 6: 1 unadjusted, 2 adj_age, 3 adj_eTIV, 4 adj_sex, 5 adj_motion, 6 fully_adjusted
#
#   USAGE
#     Rscript TWO_fda_group_difference_runs.R --list      # write/print the index table, run nothing
#     Rscript TWO_fda_group_difference_runs.R <index>     # run one model (called by the .lsf)
#
#   OUTPUT
#     results/group_differences/{metric_folder}/{config_folder}/
#
################################################################################

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

DATA_PATH    <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
RESULTS_BASE <- file.path(repo_root, "results/group_differences")
# Test overrides (leave unset for the real runs):
#   FDA_RESULTS_BASE=/some/other/folder  -> write somewhere else
#   FDA_N_BOOT=20                        -> fewer bootstrap iterations
if (nzchar(Sys.getenv("FDA_RESULTS_BASE"))) RESULTS_BASE <- Sys.getenv("FDA_RESULTS_BASE")
N_BOOT <- as.integer(Sys.getenv("FDA_N_BOOT", "1000"))
PIPELINE_SCRIPT <- file.path(repo_root, "code/03_statistical_analysis/fda_pipeline.R")

library(mgcv)   # so scat() / gaussian() resolve in the run list

# metric folder, column prefix, label, family; SW last
metrics <- list(
  list(folder = "strength", prefix = "str",               label = "Strength",           family = "gaussian"),
  list(folder = "GE_norm",  prefix = "rand_norm_wei_GE",  label = "Normalized GE",      family = "gaussian"),
  list(folder = "ACC_norm", prefix = "rand_norm_wei_ACC", label = "Normalized ACC",     family = "gaussian"),
  list(folder = "GE_raw",   prefix = "GE",                label = "Raw GE",             family = "gaussian"),
  list(folder = "ACC_raw",  prefix = "ACC",               label = "Raw ACC",            family = "gaussian"),
  list(folder = "SW",       prefix = "rand_norm_wei_SW",  label = "Small-worldness",    family = "scat")
)

# covariate configurations
configs <- list(
  list(folder = "1_unadjusted",     label = "Unadjusted",         covs = c()),
  list(folder = "2_adj_age",        label = "Adjusted: age",      covs = c("age_at_5y_mri")),
  list(folder = "3_adj_eTIV",       label = "Adjusted: eTIV",     covs = c("eTIV")),
  list(folder = "4_adj_sex",        label = "Adjusted: sex",      covs = c("sex")),
  list(folder = "5_adj_motion",     label = "Adjusted: motion",   covs = c("Rel_Motion")),
  list(folder = "6_fully_adjusted", label = "Fully adjusted",     covs = c("age_at_5y_mri", "eTIV", "sex", "Rel_Motion"))
)

# metric-specific exclusions
exclusions <- list(str = c(), rand_norm_wei_GE = c(713, 659), rand_norm_wei_ACC = c(),
                   rand_norm_wei_SW = c(), GE = c(), ACC = c())

predictor_labels <- c(Group = "Group (VPT vs FT)", age_at_5y_mri = "Age at MRI", eTIV = "eTIV", sex = "Sex",
                      Rel_Motion = "Motion (eddy relative RMS)")

# ------------------------------------------------------------------ build runs
runs <- list()
for (m in metrics) for (cf in configs) {
  preds <- c("Group", cf$covs)
  runs[[length(runs) + 1]] <- list(
    pipeline = "repol", metric = m$prefix, metric_folder = m$folder, metric_label = m$label,
    family = m$family, config_folder = cf$folder, config_label = cf$label,
    predictors = preds, skip_scaling = intersect(preds, c("Group", "sex")),
    exclude = exclusions[[m$prefix]],
    data_path = DATA_PATH,
    output_dir = file.path(RESULTS_BASE, m$folder, cf$folder))
}
stopifnot(length(runs) == 36)

index_table <- do.call(rbind, lapply(seq_along(runs), function(i) {
  r <- runs[[i]]
  data.frame(index = i, pipeline = r$pipeline, metric_folder = r$metric_folder, metric = r$metric,
             family = r$family, model = r$config_folder, predictors = paste(r$predictors, collapse = " + "),
             output_dir = r$output_dir, stringsAsFactors = FALSE)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript TWO_fda_group_difference_runs.R <index 1-36> | --list")

if (args[1] == "--list") {
  dir.create(RESULTS_BASE, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(RESULTS_BASE, "run_index_group_differences.csv")
  write.csv(index_table, f, row.names = FALSE)
  print(index_table[, c("index", "pipeline", "metric_folder", "family", "model", "predictors")], row.names = FALSE)
  cat(sprintf("\nWrote %s\n", f))
  quit(status = 0)
}

run_index <- as.integer(args[1])
if (is.na(run_index) || run_index < 1 || run_index > length(runs))
  stop(sprintf("Index %s out of range 1-%d", args[1], length(runs)))
r <- runs[[run_index]]

cat(sprintf("FDA group differences - index %d / %d\n", run_index, length(runs)))
cat(sprintf("Started %s on %s\n", format(Sys.time()), Sys.info()[["nodename"]]))
if (!file.exists(r$data_path)) stop("Data file not found: ", r$data_path)

BATCH_MODE                 <<- TRUE
BATCH_metric               <<- r$metric
BATCH_all_predictors       <<- r$predictors
BATCH_vars_to_skip_scaling <<- r$skip_scaling
BATCH_density_min          <<- 11
BATCH_density_max          <<- 100
BATCH_n_bootstrap          <<- N_BOOT
BATCH_p_threshold          <<- 0.001
BATCH_pffr_k_basis         <<- 20
BATCH_pffr_family          <<- if (r$family == "scat") scat() else gaussian()
BATCH_subjects_to_exclude  <<- r$exclude
BATCH_apply_gba_sqrt       <<- FALSE
BATCH_outlier_percentile   <<- 0.05
BATCH_output_dir           <<- r$output_dir
BATCH_data_path            <<- r$data_path
BATCH_predictor_labels     <<- predictor_labels
BATCH_run_meta             <<- list(index = run_index, pipeline = r$pipeline, analysis = "group_differences",
                                    metric_label = r$metric_label, model_label = r$config_label)

tryCatch({
  source(PIPELINE_SCRIPT)
  cat("\n=== RUN COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat(sprintf("\n!!! RUN FAILED: %s\n", conditionMessage(e)))
  quit(status = 1)
})
cat(sprintf("Finished %s\n", format(Sys.time())))