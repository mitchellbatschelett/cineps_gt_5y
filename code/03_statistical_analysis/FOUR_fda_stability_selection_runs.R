################################################################################
#
#   FOUR_fda_stability_selection_runs.R  -  stability selection
#
#   6 runs, each fit by FOUR_fda_stability_selection_pipeline.R:
#     str, rand_norm_wei_GE, rand_norm_wei_ACC, rand_norm_wei_SW (scat final family),
#     rand_norm_wei_ACC high_gba_rem, rand_norm_wei_ACC gba_binary (Supp Note 6)
#   exposures: bpd2, bw_z_new, ga, globalbrainscore2 (sqrt), anyrop, sepsis2,
#              anyivh, hydrocephalus_dc
#   forced covariates: eTIV, sex, sriskscore, age_at_5y_mri, Rel_Motion
#   metric-specific exclusions: rand_norm_wei_GE c(713, 659); none for other metrics
#
#   INDEX TABLE
#     1 ACC_norm | 2 ACC high_gba_rem | 3 ACC gba_binary | 4 strength | 5 GE_norm | 6 SW
#
#   USAGE (Terminal)
#     Rscript FOUR_fda_stability_selection_runs.R --list     # print the index table
#     Rscript FOUR_fda_stability_selection_runs.R 1          # run one
#   Run from Terminal rather than RStudio: the bootstrap forks worker processes.
#
#   OUTPUT
#     results/stability_selection/{metric_folder}/
#     results/stability_selection/ACC_sensitivity/{branch}/
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

DATA_PATH       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
OUT_BASE        <- file.path(repo_root, "results/stability_selection")
PIPELINE_SCRIPT <- file.path(repo_root, "code/03_statistical_analysis/FOUR_fda_stability_selection_pipeline.R")

library(mgcv)   # so scat() / gaussian() resolve

base_runs <- list(
  list(metric = "rand_norm_wei_ACC", folder = "ACC_norm", branch = "main",         family = "gaussian", gbs2 = "sqrt"),
  list(metric = "rand_norm_wei_ACC", folder = "ACC_norm", branch = "high_gba_rem", family = "gaussian", gbs2 = "sqrt"),
  list(metric = "rand_norm_wei_ACC", folder = "ACC_norm", branch = "gba_binary",   family = "gaussian", gbs2 = "none"),
  list(metric = "str",               folder = "strength", branch = "main",         family = "gaussian", gbs2 = "sqrt"),
  list(metric = "rand_norm_wei_GE",  folder = "GE_norm",  branch = "main",         family = "gaussian", gbs2 = "sqrt"),
  list(metric = "rand_norm_wei_SW",  folder = "SW",       branch = "main",         family = "scat",     gbs2 = "sqrt")
)

exposures_default             <- c("bpd2", "bw_z_new", "ga", "globalbrainscore2", "anyrop", "sepsis2", "anyivh", "hydrocephalus_dc")
categorical_exposures_default <- c("bpd2", "anyrop", "sepsis2", "anyivh", "hydrocephalus_dc")
categorical_forced_default    <- c("sex")

runs <- list()
for (b in base_runs) {
  out_dir <- if (b$branch == "main") {
    file.path(OUT_BASE, b$folder)
  } else {
    file.path(OUT_BASE, "ACC_sensitivity", b$branch)
  }
  runs[[length(runs) + 1]] <- c(b, list(
    pipeline  = "repol",
    forced    = c("eTIV", "sex", "sriskscore", "age_at_5y_mri", "Rel_Motion"),
    data_path = DATA_PATH,
    out_dir   = out_dir))
}
stopifnot(length(runs) == 6)

index_table <- do.call(rbind, lapply(seq_along(runs), function(i) {
  r <- runs[[i]]
  data.frame(index = i, pipeline = r$pipeline, metric = r$metric, branch = r$branch,
             final_family = r$family, forced = paste(r$forced, collapse = " + "),
             output_dir = r$out_dir, stringsAsFactors = FALSE)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript FOUR_fda_stability_selection_runs.R <index 1-6> | --list")
if (args[1] == "--list") {
  print(index_table[, c("index", "pipeline", "metric", "branch", "final_family")], row.names = FALSE)
  dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)
  write.csv(index_table, file.path(OUT_BASE, "run_index_stability_selection.csv"), row.names = FALSE)
  quit(status = 0)
}

run_index <- as.integer(args[1])
if (is.na(run_index) || run_index < 1 || run_index > length(runs))
  stop(sprintf("Index %s out of range 1-%d", args[1], length(runs)))
r <- runs[[run_index]]
if (!file.exists(r$data_path))       stop("Data file not found: ", r$data_path)
if (!file.exists(PIPELINE_SCRIPT))   stop("Pipeline not found: ", PIPELINE_SCRIPT)
dir.create(r$out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Stability selection - index %d / %d | %s | %s | %s | final family %s\n",
            run_index, length(runs), r$pipeline, r$metric, r$branch, r$family))
cat(sprintf("Started %s\n", format(Sys.time())))

BATCH_MODE                  <<- TRUE
BATCH_metric                <<- r$metric
BATCH_forced_covariates     <<- r$forced
BATCH_exposures             <<- exposures_default
BATCH_categorical_exposures <<- categorical_exposures_default
BATCH_categorical_forced    <<- categorical_forced_default
BATCH_density_min           <<- 11
BATCH_density_max           <<- 100
BATCH_pffr_family_final     <<- if (r$family == "scat") scat() else gaussian()
BATCH_subjects_to_exclude   <<- if (r$metric == "rand_norm_wei_GE") c(713, 659) else c()
BATCH_sensitivity_branch    <<- r$branch
BATCH_gbs2_transform        <<- r$gbs2
BATCH_output_dir            <<- r$out_dir
BATCH_data_path             <<- r$data_path

tryCatch({
  source(PIPELINE_SCRIPT)
  writeLines(sprintf("Finished %s", format(Sys.time())), file.path(r$out_dir, "RUN_COMPLETE.txt"))
  cat("\n=== RUN COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat(sprintf("\n!!! RUN FAILED: %s\n", conditionMessage(e)))
  quit(status = 1)
})
cat(sprintf("Finished %s\n", format(Sys.time())))
