################################################################################
#
#   FDA UNIVARIATE RUN TABLE
#
#   PURPOSE
#     Driver for THREE_fda_univariate_pipeline.R. Defines the 28 (exposure,
#     metric) cells of the fully-adjusted univariate grid reported in
#     Methods 2.7.4, Supp Note 5, Supp Table 6, and Supp Figs 11-14.
#
#   RUN TABLE (28 runs = 7 exposures x 4 primary metrics)
#     Exposures (7): bpd2, bw_z, ga, globalbrainscore2 (sqrt-transformed),
#                    anyrop, sepsis2, dwma_percent
#     Metrics (4):   str, rand_norm_wei_GE, rand_norm_wei_ACC, rand_norm_wei_SW
#     Family:        gaussian for all 28 (per Methods 2.7.4 SW caveat)
#     Forced covariates: eTIV, sex, sriskscore, age_at_5y_mri, Rel_Motion
#
#   USAGE
#     Rscript THREE_fda_univariate_runs.R <run_index>
#     where run_index is 1..28. Use submit_FDA_array_univariate.lsf to
#     submit the full set as an LSF array job.
#
#   COMPUTATIONAL COST
#     See header of THREE_fda_univariate_pipeline.R. Roughly 1-6 hours per
#     cell on a 24-core HPC node including the 1000-iteration bootstrap.
#
################################################################################

# ==============================================================================
# COMMAND-LINE INDEX
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript THREE_fda_univariate_runs.R <run_index>\n",
       "  e.g., Rscript THREE_fda_univariate_runs.R 1")
}
run_index <- as.integer(args[1])

cat(sprintf("FDA Univariate Run - Index: %d\n", run_index))
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat(sprintf("Node:      %s\n\n", Sys.info()["nodename"]))

library(mgcv)

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

pipeline_script <- file.path(repo_root, "code/03_statistical_analysis/THREE_fda_univariate_pipeline.R")
output_dir_base <- file.path(repo_root, "results/fda_univariate")
data_path       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
# --- Cores: pick ONE line below (comment out the other) ---
# n_cores <- 24                                   # fixed count (e.g. cluster node)
n_cores <- max(1, parallel::detectCores() - 1)    # auto-detect (laptop-safe default)

# ==============================================================================
# SHARED PARAMETERS
# ==============================================================================

forced_covariates_default <- c("eTIV", "sex", "sriskscore",
                               "age_at_5y_mri", "Rel_Motion")

# Binary/categorical variables that should NOT be z-scored. Everything
# else (continuous exposures of interest AND continuous forced covariates:
# eTIV, sriskscore, age_at_5y_mri, Rel_Motion) gets z-scored in the pipeline.
binary_vars <- c("sex", "bpd2", "anyrop", "sepsis2")

# Metric-specific outliers (Methods 2.7.2)
metric_outliers <- list(
  "str"               = c(128, 713),
  "rand_norm_wei_GE"  = c(659),
  "rand_norm_wei_ACC" = c(),
  "rand_norm_wei_SW"  = c(309, 321, 8155)
)

# ==============================================================================
# EXPOSURE AND METRIC LISTS
# ==============================================================================

# Order here defines exposure ordering in the run table. 7 candidate
# neonatal exposures (Methods 2.7.4).
exposure_specs <- list(
  list(name = "bpd2",              label = "BPD",    is_binary = TRUE,  needs_gba_sqrt = FALSE),
  list(name = "bw_z",              label = "BWZ",    is_binary = FALSE, needs_gba_sqrt = FALSE),
  list(name = "ga",                label = "GA",     is_binary = FALSE, needs_gba_sqrt = FALSE),
  list(name = "globalbrainscore2", label = "GBA",    is_binary = FALSE, needs_gba_sqrt = TRUE),
  list(name = "anyrop",            label = "ROP",    is_binary = TRUE,  needs_gba_sqrt = FALSE),
  list(name = "sepsis2",           label = "Sepsis", is_binary = TRUE,  needs_gba_sqrt = FALSE),
  list(name = "dwma_percent",      label = "DWMA",   is_binary = FALSE, needs_gba_sqrt = FALSE)
)

# Order here defines metric ordering in the run table. 4 primary metrics.
metric_specs <- list(
  list(metric = "str"),
  list(metric = "rand_norm_wei_GE"),
  list(metric = "rand_norm_wei_ACC"),
  list(metric = "rand_norm_wei_SW")
)

# ==============================================================================
# BUILD RUN TABLE
# ==============================================================================

# Order: outer loop = metric, inner loop = exposure
# So indices 1-7 = strength x (7 exposures), 8-14 = norm GE x (7 exposures), etc.
# This groups runs by metric, which matches the by-outcome output organization.
runs <- list()
for (m in metric_specs) {
  for (e in exposure_specs) {
    # Skip scaling for binary/categorical variables only. All continuous
    # predictors (exposures of interest AND forced covariates) get z-scored.
    skip_set <- binary_vars
    if (e$is_binary) skip_set <- unique(c(skip_set, e$name))
    
    runs[[length(runs) + 1]] <- list(
      metric                = m$metric,
      exposure              = e$name,
      forced_covariates     = forced_covariates_default,
      vars_to_skip_scaling  = skip_set,
      subjects_to_exclude   = metric_outliers[[m$metric]],
      pffr_family           = gaussian(),
      apply_gba_sqrt        = e$needs_gba_sqrt,
      density_min           = 11,
      density_max           = 100,
      n_bootstrap           = 1000,
      pffr_k_basis          = 20,
      p_threshold           = 0.001
    )
  }
}
stopifnot(length(runs) == 28)

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

# Grouped by outcome metric:
#   results/fda_univariate/{metric}/{metric}_univariate_{exposure}_11-100/
density_str    <- sprintf("%.0f-%.0f", run$density_min, run$density_max)
output_dirname <- paste(run$metric, "univariate", run$exposure,
                        density_str, sep = "_")
output_dir     <- file.path(output_dir_base, run$metric, output_dirname)

cat(sprintf("Metric:     %s\n", run$metric))
cat(sprintf("Exposure:   %s\n", run$exposure))
cat(sprintf("Outliers:   %s\n",
            if (length(run$subjects_to_exclude) > 0)
              paste(run$subjects_to_exclude, collapse = ", ") else "none"))
cat(sprintf("Output:     %s\n\n", output_dir))

# ==============================================================================
# SET BATCH VARIABLES AND SOURCE THE PIPELINE
# ==============================================================================

BATCH_MODE                 <<- TRUE
BATCH_metric               <<- run$metric
BATCH_exposure             <<- run$exposure
BATCH_forced_covariates    <<- run$forced_covariates
BATCH_vars_to_skip_scaling <<- run$vars_to_skip_scaling
BATCH_density_min          <<- run$density_min
BATCH_density_max          <<- run$density_max
BATCH_n_bootstrap          <<- run$n_bootstrap
BATCH_p_threshold          <<- run$p_threshold
BATCH_pffr_k_basis         <<- run$pffr_k_basis
BATCH_pffr_family          <<- run$pffr_family
BATCH_subjects_to_exclude  <<- run$subjects_to_exclude
BATCH_apply_gba_sqrt       <<- run$apply_gba_sqrt
BATCH_output_dir           <<- output_dir
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