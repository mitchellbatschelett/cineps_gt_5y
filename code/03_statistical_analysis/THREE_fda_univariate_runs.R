################################################################################
#
#   THREE_fda_univariate_runs.R  -  within-VPT univariate exposure models
#
#   37 runs, VPT only (N = 171), each fit by fda_pipeline.R
#     model: metric(d) ~ exposure + eTIV + sex + sriskscore + age_at_5y_mri + Rel_Motion
#     metrics: strength, normalized GE, normalized ACC, SW
#     exposures: bpd2, bw_z_new, ga, globalbrainscore2 (sqrt), anyrop, sepsis2, dwma_percent,
#                anyivh, hydrocephalus_dc
#     family: gaussian() for ALL metrics, including SW
#     z-scored: every continuous predictor (continuous exposures + eTIV, sriskscore,
#               age_at_5y_mri, Rel_Motion); not scaled: sex, bpd2, anyrop, sepsis2,
#               anyivh, hydrocephalus_dc
#     metric-specific exclusions: rand_norm_wei_GE c(713, 659); none for other metrics
#
#   INDEX TABLE
#     1-28:  1-7 strength | 8-14 GE_norm | 15-21 ACC_norm | 22-28 SW
#            within each block of 7: BPD, BWZ, GA, GBA, ROP, Sepsis, DWMA
#     29:    Supplementary Note 6 sensitivity: ACC ~ GBA with the four high-GBA
#            participants (GBA > 14) removed
#     30-37: 30-31 strength | 32-33 GE_norm | 34-35 ACC_norm | 36-37 SW
#            within each block of 2: IVH, Hydrocephalus
#
#   USAGE
#     Rscript THREE_fda_univariate_runs.R --list     # write/print the index table, run nothing
#     Rscript THREE_fda_univariate_runs.R <index>    # run one model (called by the .lsf)
#
#   OUTPUT
#     results/univariate/{metric_folder}/{exposure_folder}/
#     results/univariate_sensitivity/ACC_norm/GBA_high_gba_rem/   (run 29)
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

DATA_PATH    <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")
RESULTS_BASE <- file.path(repo_root, "results")
# Test overrides (leave unset for the real runs):
#   FDA_RESULTS_BASE=/some/other/folder  -> write somewhere else
#   FDA_N_BOOT=20                        -> fewer bootstrap iterations
if (nzchar(Sys.getenv("FDA_RESULTS_BASE"))) RESULTS_BASE <- Sys.getenv("FDA_RESULTS_BASE")
N_BOOT <- as.integer(Sys.getenv("FDA_N_BOOT", "1000"))
PIPELINE_SCRIPT <- file.path(repo_root, "code/03_statistical_analysis/fda_pipeline.R")

library(mgcv)

metrics <- list(
  list(folder = "strength", prefix = "str",               label = "Strength"),
  list(folder = "GE_norm",  prefix = "rand_norm_wei_GE",  label = "Normalized GE"),
  list(folder = "ACC_norm", prefix = "rand_norm_wei_ACC", label = "Normalized ACC"),
  list(folder = "SW",       prefix = "rand_norm_wei_SW",  label = "Small-worldness")
)

exposures <- list(
  list(folder = "BPD",    var = "bpd2",              label = "BPD"),
  list(folder = "BWZ",    var = "bw_z_new",              label = "Birth weight z-score"),
  list(folder = "GA",     var = "ga",                label = "Gestational age"),
  list(folder = "GBA",    var = "globalbrainscore2", label = "GBA (sqrt)"),
  list(folder = "ROP",    var = "anyrop",            label = "ROP"),
  list(folder = "Sepsis", var = "sepsis2",           label = "Sepsis"),
  list(folder = "DWMA",   var = "dwma_percent",      label = "DWMA (%)")
)

# cranial ultrasound exposures (runs 30-37)
exposures_cus <- list(
  list(folder = "IVH",           var = "anyivh",           label = "Any IVH"),
  list(folder = "Hydrocephalus", var = "hydrocephalus_dc", label = "Hydrocephalus")
)

binary_vars <- c("sex", "bpd2", "anyrop", "sepsis2", "anyivh", "hydrocephalus_dc")

# metric-specific exclusions
exclusions <- list(str = c(), rand_norm_wei_GE = c(713, 659), rand_norm_wei_ACC = c(), rand_norm_wei_SW = c())

predictor_labels <- c(eTIV = "eTIV", sex = "Sex", sriskscore = "Social risk score", age_at_5y_mri = "Age at MRI",
                      Rel_Motion = "Motion (eddy relative RMS)",
                      setNames(sapply(c(exposures, exposures_cus), `[[`, "label"),
                               sapply(c(exposures, exposures_cus), `[[`, "var")))

# ------------------------------------------------------------------ build runs
forced <- c("eTIV", "sex", "sriskscore", "age_at_5y_mri", "Rel_Motion")
make_run <- function(m, ex) {
  preds <- c(ex$var, forced)
  list(
    pipeline = "repol", metric = m$prefix, metric_folder = m$folder, metric_label = m$label,
    exposure_folder = ex$folder, exposure = ex$var, exposure_label = ex$label,
    predictors = preds, skip_scaling = intersect(preds, binary_vars),
    exclude = exclusions[[m$prefix]],
    data_path = DATA_PATH,
    output_dir = file.path(RESULTS_BASE, "univariate", m$folder, ex$folder))
}
runs <- list()
for (m in metrics) for (ex in exposures) runs[[length(runs) + 1]] <- make_run(m, ex)
# Supp Note 6: ACC ~ GBA with the four high-GBA participants (GBA > 14) removed
runs[[29]] <- runs[[18]]
runs[[29]]$exclude <- c(354, 559, 332, 650)
runs[[29]]$output_dir <- file.path(RESULTS_BASE, "univariate_sensitivity", "ACC_norm", "GBA_high_gba_rem")
# cranial ultrasound exposures
for (m in metrics) for (ex in exposures_cus) runs[[length(runs) + 1]] <- make_run(m, ex)
stopifnot(length(runs) == 37)

index_table <- do.call(rbind, lapply(seq_along(runs), function(i) {
  r <- runs[[i]]
  data.frame(index = i, pipeline = r$pipeline, metric_folder = r$metric_folder, metric = r$metric,
             family = "gaussian", exposure = r$exposure_folder, predictors = paste(r$predictors, collapse = " + "),
             output_dir = r$output_dir, stringsAsFactors = FALSE)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript THREE_fda_univariate_runs.R <index 1-37> | --list")

if (args[1] == "--list") {
  f <- file.path(RESULTS_BASE, "univariate", "run_index_univariate.csv")
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  write.csv(index_table, f, row.names = FALSE)
  print(index_table[, c("index", "pipeline", "metric_folder", "family", "exposure", "predictors")], row.names = FALSE)
  cat(sprintf("\nWrote %s\n", f))
  quit(status = 0)
}

run_index <- as.integer(args[1])
if (is.na(run_index) || run_index < 1 || run_index > length(runs))
  stop(sprintf("Index %s out of range 1-%d", args[1], length(runs)))
r <- runs[[run_index]]

cat(sprintf("FDA univariate - index %d / %d\n", run_index, length(runs)))
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
BATCH_pffr_family          <<- gaussian()      # gaussian for all univariate models, including SW
BATCH_subjects_to_exclude  <<- r$exclude
BATCH_apply_gba_sqrt       <<- TRUE            # only acts when globalbrainscore2 is a predictor
BATCH_outlier_percentile   <<- 0.05
BATCH_output_dir           <<- r$output_dir
BATCH_data_path            <<- r$data_path
BATCH_predictor_labels     <<- predictor_labels
BATCH_run_meta             <<- list(index = run_index, pipeline = r$pipeline, analysis = "univariate",
                                    metric_label = r$metric_label,
                                    model_label = sprintf("Exposure: %s (fully adjusted)", r$exposure_label))

tryCatch({
  source(PIPELINE_SCRIPT)
  cat("\n=== RUN COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat(sprintf("\n!!! RUN FAILED: %s\n", conditionMessage(e)))
  quit(status = 1)
})
cat(sprintf("Finished %s\n", format(Sys.time())))