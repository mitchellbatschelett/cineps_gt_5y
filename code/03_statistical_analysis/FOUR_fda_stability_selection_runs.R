################################################################################
#
#   FDA STABILITY SELECTION RUN TABLE
#
#   PURPOSE
#     Driver for FOUR_fda_stability_selection_pipeline.R. Defines the 6 runs
#     reported in the manuscript and supplement: four main-text metric runs
#     (strength, normalized GE/ACC/SW) and two ACC sensitivity branches
#     (Supp Note 6 / Supp Figs 15-16). Sets BATCH_* variables for one run
#     selected by index and sources the pipeline.
#
#   RUN TABLE (6 runs)
#     1. str                main           gaussian
#     2. rand_norm_wei_GE   main           gaussian
#     3. rand_norm_wei_ACC  main           gaussian
#     4. rand_norm_wei_SW   main           scat
#     5. rand_norm_wei_ACC  high_gba_rem   gaussian   (Supp Fig 15)
#     6. rand_norm_wei_ACC  gba_binary     gaussian   (Supp Fig 16)
#
#   USAGE
#     Rscript FOUR_fda_stability_selection_runs.R <run_index>
#     where run_index is 1..6. To run interactively, set run_index manually
#     and source.
#
#   COMPUTATIONAL COST
#     See header of FOUR_fda_stability_selection_pipeline.R. Roughly two days 
#     of total wall time for all six runs on a recent workstation.
#
################################################################################

# ==============================================================================
# COMMAND-LINE INDEX
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript FOUR_fda_stability_selection_runs.R <run_index>\n",
       "  e.g., Rscript FOUR_fda_stability_selection_runs.R 1")
}
run_index <- as.integer(args[1])

cat(sprintf("FDA Stability Selection Run - Index: %d\n", run_index))
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat(sprintf("Node:      %s\n\n", Sys.info()["nodename"]))

library(mgcv)  # for gaussian() / scat() inside the runs list

# ==============================================================================
# PATHS
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

pipeline_script <- file.path(repo_root,
                             "code/03_statistical_analysis",
                             "FOUR_fda_stability_selection_pipeline.R")
data_path       <- file.path(repo_root,
                             "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
output_dir_base <- file.path(repo_root,
                             "results/fda_stability_selection")

# ==============================================================================
# SHARED PARAMETERS
# ==============================================================================

forced_covariates_default <- c("eTIV", "sex", "sriskscore",
                               "age_at_5y_mri", "Rel_Motion")
exposures_default <- c("bpd2", "bw_z", "ga", "globalbrainscore2",
                       "anyrop", "sepsis2")
categorical_exposures_default <- c("bpd2", "anyrop", "sepsis2")
categorical_forced_default    <- c("sex")

# Methods 2.7.2 metric-specific outlier exclusions
metric_outliers <- list(
  "str"               = c(128, 713),
  "rand_norm_wei_GE"  = c(659),
  "rand_norm_wei_ACC" = c(),
  "rand_norm_wei_SW"  = c(309, 321, 8155)
)

# ==============================================================================
# RUN TABLE
# ==============================================================================

runs <- list(
  
  # 1. Strength (main)
  list(metric              = "str",
       sensitivity_branch  = "main",
       pffr_family_final   = gaussian(),
       subjects_to_exclude = metric_outliers[["str"]],
       gbs2_transform      = "sqrt"),
  
  # 2. Normalized GE (main)
  list(metric              = "rand_norm_wei_GE",
       sensitivity_branch  = "main",
       pffr_family_final   = gaussian(),
       subjects_to_exclude = metric_outliers[["rand_norm_wei_GE"]],
       gbs2_transform      = "sqrt"),
  
  # 3. Normalized ACC (main; Fig 3)
  list(metric              = "rand_norm_wei_ACC",
       sensitivity_branch  = "main",
       pffr_family_final   = gaussian(),
       subjects_to_exclude = metric_outliers[["rand_norm_wei_ACC"]],
       gbs2_transform      = "sqrt"),
  
  # 4. Normalized SW (main; Fig 4) - scat family per Methods 2.7.3/2.7.4
  list(metric              = "rand_norm_wei_SW",
       sensitivity_branch  = "main",
       pffr_family_final   = scat(),
       subjects_to_exclude = metric_outliers[["rand_norm_wei_SW"]],
       gbs2_transform      = "sqrt"),
  
  # 5. ACC sensitivity: high-GBA participants removed (Supp Fig 15)
  list(metric              = "rand_norm_wei_ACC",
       sensitivity_branch  = "high_gba_rem",
       pffr_family_final   = gaussian(),
       subjects_to_exclude = metric_outliers[["rand_norm_wei_ACC"]],
       gbs2_transform      = "sqrt"),
  
  # 6. ACC sensitivity: GBA recoded as binary (>=8 vs <8) (Supp Fig 16)
  list(metric              = "rand_norm_wei_ACC",
       sensitivity_branch  = "gba_binary",
       pffr_family_final   = gaussian(),
       subjects_to_exclude = metric_outliers[["rand_norm_wei_ACC"]],
       gbs2_transform      = "none")
)
stopifnot(length(runs) == 6)

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

density_str <- "11-100"
if (run$sensitivity_branch == "main") {
  output_dirname <- paste(run$metric, "stabsel", density_str, sep = "_")
  output_dir <- file.path(output_dir_base, output_dirname)
} else {
  output_dirname <- paste(run$metric, "stabsel", density_str,
                          run$sensitivity_branch, sep = "_")
  output_dir <- file.path(output_dir_base, "ACC_sensitivity", output_dirname)
}

cat(sprintf("Metric:     %s\n", run$metric))
cat(sprintf("Sensitivity: %s\n", run$sensitivity_branch))
cat(sprintf("Family:     %s\n",
            if (inherits(run$pffr_family_final, "general.family"))
              "scat" else run$pffr_family_final$family))
cat(sprintf("Outliers:   %s\n",
            if (length(run$subjects_to_exclude) > 0)
              paste(run$subjects_to_exclude, collapse = ", ") else "none"))
cat(sprintf("Output:     %s\n\n", output_dir))

# ==============================================================================
# SET BATCH VARIABLES AND SOURCE THE PIPELINE
# ==============================================================================

BATCH_MODE                  <<- TRUE
BATCH_metric                <<- run$metric
BATCH_forced_covariates     <<- forced_covariates_default
BATCH_exposures             <<- exposures_default
BATCH_categorical_exposures <<- categorical_exposures_default
BATCH_categorical_forced    <<- categorical_forced_default
BATCH_density_min           <<- 11
BATCH_density_max           <<- 100
BATCH_pffr_family_final     <<- run$pffr_family_final
BATCH_subjects_to_exclude   <<- run$subjects_to_exclude
BATCH_sensitivity_branch    <<- run$sensitivity_branch
BATCH_gbs2_transform        <<- run$gbs2_transform
BATCH_output_dir            <<- output_dir
BATCH_data_path             <<- data_path

tryCatch({
  source(pipeline_script)
  cat("\n\n=== RUN COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat(sprintf("\n\n!!! RUN FAILED: %s\n", e$message))
  quit(status = 1)
})

cat(sprintf("Finished: %s\n", Sys.time()))