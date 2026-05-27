# ============================================================
# TWO_merge_clinical_and_metrics.R
#
# PURPOSE: Merge per-subject graph theory metrics with clinical and
#          demographic data to produce the analysis-ready files used by
#          all downstream statistical analyses.
#
#          Produces two outputs:
#            1. cohort_171VPT_45FT_postVQC.xlsx -- all subjects, Group included
#            2. cohort_171VPT_postVQC.xlsx     -- VPT subjects only, Group dropped
#
# INPUTS:
#   - GT_PATH:       data/intermediate/graph_theory_metrics_per_subject.xlsx
#                    (output of ONE_compute_graph_theory_metrics.m)
#   - CLINICAL_PATH: data/demographic_clinical/demo_clinical_only.xlsx
#
# OUTPUTS:
#   - OUT_FULL:      data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx
#   - OUT_VPT:       data/analysis_ready/cohort_171VPT_postVQC.xlsx
#
# REQUIRES: R packages readxl, writexl
# ============================================================

# ============================================================
# CONFIG -- edit before running
# ============================================================
# Repo root is located automatically so a fresh clone runs without edits.
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

# GT_PATH is the per-subject metrics workbook written by the 02_graph_theory
# step ONE script (graph_theory_metrics_per_subject.xlsx). It is provided in
# data/intermediate/; regenerate it by running ONE_compute_graph_theory_metrics.m
# if desired. Alternatively, skip this script entirely and use the already-merged
# cohort files in data/analysis_ready/.
GT_PATH       <- file.path(repo_root, "data/intermediate/graph_theory_metrics_per_subject.xlsx")
CLINICAL_PATH <- file.path(repo_root, "data/demographic_clinical/demo_clinical_only.xlsx")
OUT_FULL      <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
OUT_VPT       <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_postVQC.xlsx")

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
})

# --- Load data ---
gt       <- read_excel(GT_PATH)
clinical <- read_excel(CLINICAL_PATH)

cat(sprintf("GT file:       %d subjects, %d columns\n", nrow(gt), ncol(gt)))
cat(sprintf("Clinical file: %d subjects, %d columns\n", nrow(clinical), ncol(clinical)))

# --- Harmonize IDs ---
# GT file (from MATLAB) writes IDs as "sub-###"; clinical uses numeric IDs.
gt$ID       <- as.numeric(gsub("^sub-", "", gt$ID))
clinical$ID <- as.numeric(clinical$ID)

# --- Sanity check: drop Group from GT if present, since clinical is source of truth ---
if ("Group" %in% names(gt)) {
  gt$Group <- NULL
}

# --- Merge on ID (inner join) ---
merged <- merge(clinical, gt, by = "ID", all = FALSE)
cat(sprintf("Merged:        %d subjects, %d columns\n\n", nrow(merged), ncol(merged)))

# --- Reconciliation report ---
gt_only   <- setdiff(gt$ID, clinical$ID)
clin_only <- setdiff(clinical$ID, gt$ID)

if (length(gt_only) > 0) {
  cat(sprintf("WARNING: %d subjects in GT file but not in clinical:\n  %s\n\n",
              length(gt_only), paste(gt_only, collapse = ", ")))
}
if (length(clin_only) > 0) {
  cat(sprintf("WARNING: %d subjects in clinical file but not in GT:\n  %s\n\n",
              length(clin_only), paste(clin_only, collapse = ", ")))
}

n_vpt <- sum(merged$Group == 1, na.rm = TRUE)
n_ft  <- sum(merged$Group == 0, na.rm = TRUE)
cat(sprintf("Merged cohort breakdown: %d VPT + %d FT = %d total\n\n",
            n_vpt, n_ft, n_vpt + n_ft))

# --- Ensure output directory exists ---
out_dir <- dirname(OUT_FULL)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Write full cohort (VPT + FT) ---
write_xlsx(merged, OUT_FULL)
cat(sprintf("Wrote: %s\n", OUT_FULL))

# --- Write VPT-only (drop Group) ---
vpt_only <- merged[merged$Group == 1, ]
vpt_only$Group <- NULL
write_xlsx(vpt_only, OUT_VPT)
cat(sprintf("Wrote: %s  [%d subjects]\n", OUT_VPT, nrow(vpt_only)))

cat("\nDone.\n")