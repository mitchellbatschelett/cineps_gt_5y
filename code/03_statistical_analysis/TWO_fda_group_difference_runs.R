################################################################################
#
#   FDA GROUP DIFFERENCE RUN TABLE
#
#   PURPOSE
#     Driver for TWO_fda_pipeline.R. Defines the 36 (metric x covariate
#     configuration) cells reported in the manuscript and supplement, sets
#     BATCH_* variables for one run selected by a command-line index, and
#     sources the pipeline.
#
#   RUN TABLE (36 runs = 6 metrics x 6 covariate configurations)
#     Metrics (6):    str, rand_norm_wei_GE, rand_norm_wei_ACC,
#                     rand_norm_wei_SW, GE (raw), ACC (raw)
#     Configurations (6):
#       1. Group                                                  [unadjusted]
#       2. Group + age_at_5y_mri
#       3. Group + eTIV
#       4. Group + sex
#       5. Group + Rel_Motion
#       6. Group + age_at_5y_mri + eTIV + sex + Rel_Motion        [fully adjusted]
#
#     Family: gaussian() for all metrics except rand_norm_wei_SW, which uses
#     scat() per Methods 2.7.3.
#
#     Metric-specific outlier exclusions (Methods 2.7.2):
#       str:               c(128, 713)
#       rand_norm_wei_GE:  c(659)
#       rand_norm_wei_ACC: c()
#       rand_norm_wei_SW:  c(309, 321, 8155)
#       GE  (raw):         c(128, 713)
#       ACC (raw):         c(128, 713)
#
#   USAGE
#     Rscript TWO_fda_group_difference_runs.R <run_index>
#     where run_index is 1..36. Use submit_FDA_array_group_diff.lsf to submit
#     the full set as an LSF array job, if using HPC cluster. To run interactively, 
#     set run_index manually and source this file.
#
#   COMPUTATIONAL COST
#     See header of TWO_fda_pipeline.R. Hours per gaussian run; much longer
#     per scat-SW run (runs 19-24).
#
################################################################################

# ==============================================================================
# COMMAND-LINE INDEX
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript TWO_fda_group_difference_runs.R <run_index>\n",
       "  e.g., Rscript TWO_fda_group_difference_runs.R 1")
}
run_index <- as.integer(args[1])

cat(sprintf("FDA Group Difference Run - Index: %d\n", run_index))
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat(sprintf("Node:      %s\n\n", Sys.info()["nodename"]))

# mgcv must be loaded so that scat() resolves before being placed into the runs list
library(mgcv)

# ==============================================================================
# PATHS - cluster defaults; uncomment LOCAL OVERRIDE block for laptop use
# ==============================================================================

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

pipeline_script <- file.path(repo_root, "code/03_statistical_analysis/TWO_fda_pipeline.R")
output_dir_base <- file.path(repo_root, "results/fda_group_differences")
data_path       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
# --- Cores: pick ONE line below (comment out the other) ---
# n_cores <- 24                                   # fixed count (e.g. cluster node)
n_cores <- max(1, parallel::detectCores() - 1)    # auto-detect (laptop-safe default)

# ==============================================================================
# METRIC SPECIFICATIONS
# ==============================================================================

# Order here defines the metric ordering in the run table (rows 1-6 are
# the first metric across all 6 configs, rows 7-12 are the second, etc.).
# Methods 2.7.2 outlier exclusions are hardcoded per metric here so the runs
# list below doesn't need to repeat them.
metric_specs <- list(
  list(metric = "str",
       family = gaussian(),
       exclude = c(128, 713)),
  list(metric = "rand_norm_wei_GE",
       family = gaussian(),
       exclude = c(659)),
  list(metric = "rand_norm_wei_ACC",
       family = gaussian(),
       exclude = c()),
  list(metric = "rand_norm_wei_SW",
       family = scat(),
       exclude = c(309, 321, 8155)),
  list(metric = "GE",
       family = gaussian(),
       exclude = c(128, 713)),
  list(metric = "ACC",
       family = gaussian(),
       exclude = c(128, 713))
)

# ==============================================================================
# COVARIATE CONFIGURATIONS
# ==============================================================================

# Order here defines covariate-config ordering within each metric block.
# vars_to_skip_scaling encodes "do not z-score before fitting".
covar_configs <- list(
  list(predictors = c("Group"),
       skip_scaling = c("Group")),
  list(predictors = c("Group", "age_at_5y_mri"),
       skip_scaling = c("Group")),
  list(predictors = c("Group", "eTIV"),
       skip_scaling = c("Group")),
  list(predictors = c("Group", "sex"),
       skip_scaling = c("Group", "sex")),
  list(predictors = c("Group", "Rel_Motion"),
       skip_scaling = c("Group")),
  list(predictors = c("Group", "age_at_5y_mri", "eTIV", "sex", "Rel_Motion"),
       skip_scaling = c("Group", "sex"))
)

# ==============================================================================
# BUILD RUN TABLE (cross metric_specs x covar_configs)
# ==============================================================================

runs <- list()
for (m in metric_specs) {
  for (c in covar_configs) {
    runs[[length(runs) + 1]] <- list(
      metric               = m$metric,
      pffr_family          = m$family,
      subjects_to_exclude  = m$exclude,
      predictors           = c$predictors,
      vars_to_skip_scaling = c$skip_scaling,
      density_min          = 11,
      density_max          = 100,
      n_bootstrap          = 1000,
      pffr_k_basis         = 20,
      p_threshold          = 0.001
    )
  }
}
stopifnot(length(runs) == 36)

# ==============================================================================
# VALIDATE INDEX
# ==============================================================================

if (run_index < 1 || run_index > length(runs)) {
  stop(sprintf("Run index %d is out of range. Valid range: 1-%d",
               run_index, length(runs)))
}

run <- runs[[run_index]]

# ==============================================================================
# BUILD OUTPUT DIRECTORY NAME
# ==============================================================================

# Convention: {metric}_FDA_{predictors_underscored}_{dmin}-{dmax}_fullsample
# (preserved from the upstream pipeline so downstream figure scripts find it)
pred_str       <- paste(run$predictors, collapse = "_")
density_str    <- sprintf("%.0f-%.0f", run$density_min, run$density_max)
output_dirname <- paste(run$metric, "FDA", pred_str, density_str,
                        "fullsample", sep = "_")

cat(sprintf("Metric:     %s\n", run$metric))
cat(sprintf("Family:     %s\n", run$pffr_family$family))
cat(sprintf("Predictors: %s\n", paste(run$predictors, collapse = ", ")))
cat(sprintf("Outliers:   %s\n",
            if (length(run$subjects_to_exclude) > 0)
              paste(run$subjects_to_exclude, collapse = ", ") else "none"))
cat(sprintf("Output:     %s\n\n", output_dirname))

# ==============================================================================
# SET BATCH VARIABLES AND SOURCE THE PIPELINE
# ==============================================================================

BATCH_MODE                 <<- TRUE
BATCH_metric               <<- run$metric
BATCH_all_predictors       <<- run$predictors
BATCH_vars_to_skip_scaling <<- run$vars_to_skip_scaling
BATCH_density_min          <<- run$density_min
BATCH_density_max          <<- run$density_max
BATCH_n_bootstrap          <<- run$n_bootstrap
BATCH_p_threshold          <<- run$p_threshold
BATCH_pffr_k_basis         <<- run$pffr_k_basis
BATCH_pffr_family          <<- run$pffr_family
BATCH_subjects_to_exclude  <<- run$subjects_to_exclude
BATCH_output_dir           <<- file.path(output_dir_base, output_dirname)
BATCH_data_path            <<- data_path
BATCH_n_cores              <<- n_cores

tryCatch({
  source(pipeline_script)
  cat("\n\n=== RUN COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat(sprintf("\n\n!!! RUN FAILED: %s\n", e$message))
  quit(status = 1)
})

cat(sprintf("Finished: %s\n", Sys.time()))