################################################################################
#                                                                              #
#   CINEPS GT MANUSCRIPT - MAIN RESULTS FIGURES & TABLES                       #
#                                                                              #
#   Generates the main-text figures (2-5) and tables (1-3) for the manuscript. #
#   Figure 1 is produced directly by ONE_density_group_difference.R and is    #
#   referenced here only for completeness; no rebuild logic is included.       #
#                                                                              #
#   STRUCTURE                                                                  #
#     Section 0:  Configuration (paths, libraries, plot theme, STRICT flag)    #
#     Section 1:  Helpers (RData loading, beta extraction, CI utilities)       #
#     Section 2:  Data loading + preprocessing                                 #
#     Section 3:  Table 1 - Participant Characteristics                        #
#     Section 4:  Table 2 + Figure 2 - VPT vs FT Group Differences             #
#     Section 5:  Table 3 + Figures 3-5 - Stability Selection                  #
#                                                                              #
#   Each section is self-contained and can be run in isolation after Sections  #
#   0-2 are sourced.                                                           #
#                                                                              #
#   DEPENDENCIES                                                               #
#     - data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx                    #
#     - results/density_group_difference/density_sensitivity_results.xlsx      #
#     - results/fda_group_differences/{metric}_FDA_{...}/                      #
#         {metric}_FDA_results.RData  (12 directories, 4 metrics x 2 configs)  #
#     - results/fda_stability_selection/{metric}_stabsel_11-100/               #
#         {metric}_stabsel_results.RData  (4 main-run directories)             #
#                                                                              #
#   OUTPUTS                                                                    #
#     - tables/main/Table_1_participant_characteristics.csv                   #
#     - tables/main/Table_2_group_differences.csv                             #
#     - tables/main/Table_3_stability_selection.csv                           #
#     - figures/main/Figure_2_group_difference_betas.{pdf,png}                 #
#     - figures/main/Figure_3_ACC_stability_selection.{pdf,png}                #
#     - figures/main/Figure_4_SW_stability_selection.{pdf,png}                 #
#     - figures/main/Figure_5_GE_strength_stability.{pdf,png}                  #
#                                                                              #
#   REQUIRES                                                                   #
#     R 4.4.0; readxl, dplyr, tidyr, tibble, ggplot2, patchwork, refund, mgcv  #
#                                                                              #
################################################################################


# ============================================================================ 
# SECTION 0: CONFIGURATION
# ============================================================================ 

# --- STRICT flag --------------------------------------------------------------
# TRUE  (ship): any missing input file (RData, xlsx) is a fatal error. This is
#               the configuration for the repository and for reproducibility
#               reviewers running the pipeline end-to-end.
# FALSE (dev):  missing inputs are reported and the affected section is skipped.

STRICT <- TRUE

# --- Repo root and derived paths ---------------------------------------------
# Locate the repository root automatically so a fresh clone runs without edits.
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

# Demographic / clinical variables only - used by Table 1. Avoids loading the
# 2800-column merged file just to compute participant characteristics.
demo_path  <- file.path(repo_root, "data/demographic_clinical/demo_clinical_only.xlsx")

# Full analysis-ready cohort (clinical + GT density columns) - used by any
# section that needs density-resolved metrics. Not needed by Table 1.
data_path  <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
density_xlsx <- file.path(repo_root,
                          "results/density_group_difference/density_sensitivity_results.xlsx")

fda_group_root <- file.path(repo_root, "results/fda_group_differences")
fda_stabsel_root <- file.path(repo_root, "results/fda_stability_selection")

figures_out <- file.path(repo_root, "figures/main")
tables_out  <- file.path(repo_root, "tables/main")

for (d in c(figures_out, tables_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Analytic density range (Methods 2.7.2) ----------------------------------
density_min <- 11
density_max <- 100

# --- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)    # read cohort xlsx + density_sensitivity_results.xlsx
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(refund)    # pffr machinery; needed so plot.pffr dispatches
  library(mgcv)      # gaussian / scat families
})

# --- Common plot theme -------------------------------------------------------
theme_main <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title    = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30"),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 2, color = "black"),
      legend.title  = element_text(size = base_size - 1, face = "bold"),
      legend.text   = element_text(size = base_size - 2)
    )
}

# --- Save utility: PNG @ 300 dpi + vector PDF --------------------------------
save_fig <- function(plot_obj, name, width_in, height_in, out_dir = figures_out) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(name, ".png")),
         plot = plot_obj, width = width_in, height = height_in, dpi = 300,
         units = "in", bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")),
         plot = plot_obj, width = width_in, height = height_in,
         units = "in", device = cairo_pdf)
  cat(sprintf("  Saved: %s.png, %s.pdf\n", name, name))
  invisible(plot_obj)
}

# --- Missing-input handler ---------------------------------------------------
# Centralizes STRICT vs dev behavior. Returns invisibly TRUE if file is present
# and FALSE (after warning) if it is missing and STRICT is FALSE; stops if
# STRICT is TRUE.
require_input <- function(path, label = NULL) {
  if (file.exists(path)) return(invisible(TRUE))
  msg <- if (is.null(label)) sprintf("Missing input: %s", path)
  else                sprintf("Missing input (%s): %s", label, path)
  if (STRICT) stop(msg, call. = FALSE)
  cat(sprintf("  [skip] %s\n", msg))
  invisible(FALSE)
}


# ============================================================================ 
# SECTION 1: HELPERS
# ============================================================================ 

# --- Format utilities --------------------------------------------------------

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_p <- Vectorize(fmt_p)

# As fmt_p but always emits the relational sign explicitly:
#   p < 0.001  ->  "<0.001"
#   p >= 0.001 ->  "0.025" (no equals sign; same digits as fmt_p)
# Kept as a separate function for clarity at call sites that build composite
# strings like "F = 11.95, p <0.001".
#
# Pass-through for already-formatted character p values: some inputs arrive as
# character strings (e.g. "<0.001") when read from an xlsx where Excel mixed
# numeric and string formatting in one column. In that case, return as-is.
fmt_p_relop <- function(p) {
  if (length(p) == 0) return(NA_character_)
  if (is.character(p)) {
    # Try to parse as number; if it parses, format normally. If not (e.g.
    # "<0.001"), return the string as-is.
    p_num <- suppressWarnings(as.numeric(p))
    if (is.na(p_num)) return(p)
    p <- p_num
  }
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_p_relop <- Vectorize(fmt_p_relop)

fmt_stat <- function(x, digits = 2) {
  if (is.na(x)) return(NA_character_)
  sprintf(paste0("%.", digits, "f"), x)
}
fmt_stat <- Vectorize(fmt_stat, vectorize.args = "x")

fmt_pct <- function(x, digits = 1) {
  if (is.na(x)) return(NA_character_)
  sprintf(paste0("%.", digits, "f%%"), x * 100)
}
fmt_pct <- Vectorize(fmt_pct, vectorize.args = "x")

# Format a variance-explained value as a percentage. Accepts either a fraction
# (in [0, 1]) or a percentage (>1, <=100). NA -> em-dash. Optional suffix is
# appended to support footnote markers (e.g. superscript 'a' for scat metrics).
fmt_ve <- function(x, suffix = "") {
  if (is.na(x)) return("\u2014")
  pct <- if (x <= 1) x * 100 else x
  paste0(sprintf("%.1f%%", pct), suffix)
}

# --- PFFR plotting helpers (getPlotObject / getCIsList) -----------------------
# Reconstruct the same beta(density) curves and bootstrap CIs the per-run
# pipeline plots produce, so the assembled panels match the individual outputs.

get_plot_object <- function(model) {
  # plot.gam() needs an open graphics device to return plot coordinates; use a
  # throwaway svg() device. The captured plot data is independent of device choice.
  ff <- tempfile()
  svg(filename = ff)
  po <- plot(model)
  dev.off()
  unlink(ff)
  po
}

get_cis_list <- function(coefboot_bs) {
  smList <- coefboot_bs$smterms
  out <- vector("list", length(smList))
  for (i in seq_along(smList)) {
    out[[i]] <- list(
      x        = smList[[i]][[2]],
      y        = smList[[i]][[1]],
      ci_lower = smList[[i]][["2.5%"]],
      ci_upper = smList[[i]][["97.5%"]]
    )
  }
  out
}

# --- Load RData into a fresh environment -------------------------------------
# Avoids polluting the global namespace and lets a caller pull only the objects
# it needs. Returns the environment.
load_rdata <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e
}

# --- Find a single smooth term by name ---------------------------------------
# Term names in bootstrap_coefs$smterms look like "Group(densities)",
# "globalbrainscore2(densities)", etc. Returns the 1-based index of the first
# matching term, or NA_integer_ if none match.
find_term_index <- function(term_names, predictor) {
  # Strip "(...)" suffix for matching
  bare <- sub("\\(.*\\)$", "", term_names)
  idx <- which(bare == predictor)
  if (length(idx) == 0) {
    # Fallback: grep prefix match (handles edge cases like factor terms)
    idx <- grep(paste0("^", predictor), term_names)
  }
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

# --- Extract beta(density) curve + bootstrap CI for one predictor ------------
# Returns a data.frame with columns: density, beta, ci_lower, ci_upper.
extract_beta_curve <- function(env, predictor,
                               model_obj_name = "pffr_fit") {
  model <- env[[model_obj_name]]
  if (is.null(model)) stop(sprintf("No %s in environment", model_obj_name))
  if (is.null(env$bootstrap_coefs))
    stop("No bootstrap_coefs in environment")
  
  term_names <- names(env$bootstrap_coefs$smterms)
  idx <- find_term_index(term_names, predictor)
  if (is.na(idx)) {
    stop(sprintf("Predictor '%s' not found among terms: %s",
                 predictor, paste(term_names, collapse = ", ")))
  }
  
  po <- get_plot_object(model)
  CIs <- get_cis_list(env$bootstrap_coefs)
  
  data.frame(
    density  = po[[idx]]$x,
    beta     = po[[idx]]$fit,
    ci_lower = CIs[[idx]]$ci_lower,
    ci_upper = CIs[[idx]]$ci_upper
  )
}

# --- Compute CI exclusion range ---------------------------------------------
# Given a data.frame with density / ci_lower / ci_upper columns, returns a
# formatted string like "45.2–100%", "11–27.2%", or "ns" (if zero is contained
# in the CI at every density). Uses an en-dash for ranges per manuscript style.
#
# Behavior:
#   - If CI excludes zero (lower > 0 OR upper < 0) at every density, returns
#     "{dmin}-{dmax}%".
#   - If CI excludes zero over one or more contiguous intervals, returns the
#     union joined by ", " (e.g. "11-27.2%, 65-78%"). In practice manuscript
#     reports a single interval per panel.
#   - If CI contains zero everywhere, returns "ns".
ci_exclusion_range <- function(beta_df) {
  if (nrow(beta_df) == 0) return(NA_character_)
  d  <- beta_df$density
  lo <- beta_df$ci_lower
  hi <- beta_df$ci_upper
  excludes <- (lo > 0) | (hi < 0)
  if (!any(excludes, na.rm = TRUE)) return("ns")
  # Identify contiguous runs of TRUE
  rle_obj <- rle(excludes)
  ends   <- cumsum(rle_obj$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  segs <- mapply(function(s, e, v) {
    if (!isTRUE(v)) return(NULL)
    sprintf("%s–%s%%", fmt_density(d[s]), fmt_density(d[e]))
  }, starts, ends, rle_obj$values, SIMPLIFY = FALSE)
  segs <- Filter(Negate(is.null), segs)
  paste(unlist(segs), collapse = ", ")
}

# Format a density value: integer if it rounds to integer, else one decimal.
fmt_density <- function(d) {
  if (abs(d - round(d)) < 1e-6) return(sprintf("%d", as.integer(round(d))))
  sprintf("%.1f", d)
}
fmt_density <- Vectorize(fmt_density, vectorize.args = "d")


# ============================================================================ 
# SECTION 2: DATA LOADING + PREPROCESSING
# ============================================================================ 
# We load only the demographic / clinical slice here. Table 1 needs nothing
# beyond these columns. Sections that need density-resolved GT metrics will
# load the full cohort xlsx (data_path) on demand.

cat("\n========================================\n")
cat("Loading demographic / clinical data\n")
cat("========================================\n\n")

require_input(demo_path, "demographic xlsx")
df <- as.data.frame(read_excel(demo_path))
cat(sprintf("Loaded: %d rows x %d cols\n", nrow(df), ncol(df)))

# Group factor (Methods 2.1: Group == 1 -> VPT, Group == 0 -> FT)
df$Group <- factor(df$Group, levels = c(0, 1), labels = c("FT", "VPT"))
cat(sprintf("Group counts: %s\n",
            paste(sprintf("%s = %d", levels(df$Group), table(df$Group)),
                  collapse = ", ")))


# ============================================================================ 
# SECTION 3: TABLE 1 - PARTICIPANT CHARACTERISTICS
# ============================================================================ 
# Manuscript Table 1 (rows in order):
#   Shared variables:
#     GA at birth, weeks                Welch's t
#     Sex (Male / Female counts)        Chi-square
#     Corrected age at MRI, years       Welch's t
#     Total Intracranial Volume, cm^3   Welch's t
#     Relative Motion, mm               Welch's t
#   VPT clinical variables (VPT only; FT column = "-"):
#     Birth weight, grams               mean +/- SD
#     Maternal social risk score        mean +/- SD
#     Global brain abnormality score    median [IQR]
#     BPD grade (No, 1, 2, 3)           count (%)
#     ROP                               count (%)
#     Sepsis                            count (%)
# ============================================================================ 

cat("\n========================================\n")
cat("Building Table 1 (Participant Characteristics)\n")
cat("========================================\n\n")

# Helper: Welch's t result formatted as desired
welch_row <- function(varname, label, units, digits = 1) {
  vpt <- df[[varname]][df$Group == "VPT"]
  ft  <- df[[varname]][df$Group == "FT"]
  vpt <- vpt[!is.na(vpt)]
  ft  <- ft[!is.na(ft)]
  if (length(vpt) == 0 || length(ft) == 0) {
    return(list(VPT = NA_character_, FT = NA_character_, p = NA_character_))
  }
  tt <- t.test(vpt, ft)
  fmt <- paste0("%.", digits, "f")
  list(
    label = sprintf("%s%s", label, if (nzchar(units)) sprintf(", %s", units) else ""),
    VPT   = sprintf(paste0(fmt, " \u00b1 ", fmt), mean(vpt), sd(vpt)),
    FT    = sprintf(paste0(fmt, " \u00b1 ", fmt), mean(ft),  sd(ft)),
    p     = fmt_p(tt$p.value)
  )
}

# Helper: median [IQR] for a numeric variable (one group only)
median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  q <- quantile(x, c(0.25, 0.5, 0.75))
  fmt <- paste0("%.", digits, "f")
  sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), q[2], q[1], q[3])
}

# Helper: count (pct) formatted "n (xx.x)"
n_pct <- function(x, total) {
  if (is.na(x) || is.na(total) || total == 0) return("-")
  sprintf("%d (%.1f)", as.integer(x), 100 * x / total)
}

n_vpt <- sum(df$Group == "VPT")
n_ft  <- sum(df$Group == "FT")

# --- Build rows -------------------------------------------------------------

t1_rows <- list()

# Shared variables header (in xlsx: italic / bold via styling)
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared",
                                       Variable = "Shared variables",
                                       VPT = "", FT = "", p = "")

# GA at birth
ga_row <- welch_row("ga", "GA at birth", "weeks", digits = 1)
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared",
                                       Variable = paste0(ga_row$label, "\u1d43"),
                                       VPT = ga_row$VPT, FT = ga_row$FT,
                                       p = ga_row$p)

# Sex (chi-square)
sex_tab <- table(df$Group, df$sex)
sex_chi <- suppressWarnings(chisq.test(sex_tab))
# Determine which sex code corresponds to Female (Methods 2.1: 1 = F, 0 = M)
n_male_vpt   <- sum(df$sex == 0 & df$Group == "VPT", na.rm = TRUE)
n_female_vpt <- sum(df$sex == 1 & df$Group == "VPT", na.rm = TRUE)
n_male_ft    <- sum(df$sex == 0 & df$Group == "FT",  na.rm = TRUE)
n_female_ft  <- sum(df$sex == 1 & df$Group == "FT",  na.rm = TRUE)

t1_rows[[length(t1_rows) + 1]] <- list(section = "shared",
                                       Variable = "Sex\u1d47",
                                       VPT = "", FT = "",
                                       p = fmt_p(sex_chi$p.value))
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared_sub",
                                       Variable = "  Male",
                                       VPT = n_pct(n_male_vpt, n_vpt),
                                       FT  = n_pct(n_male_ft, n_ft),
                                       p = "")
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared_sub",
                                       Variable = "  Female",
                                       VPT = n_pct(n_female_vpt, n_vpt),
                                       FT  = n_pct(n_female_ft, n_ft),
                                       p = "")

# Corrected age at MRI
age_row <- welch_row("age_at_5y_mri", "Corrected age at MRI", "years", digits = 2)
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared",
                                       Variable = paste0(age_row$label, "\u1d43"),
                                       VPT = age_row$VPT, FT = age_row$FT,
                                       p = age_row$p)

# TIV — eTIV in mm^3, convert to cm^3 (divide by 1000) per manuscript
tiv_vpt <- df$eTIV[df$Group == "VPT"] / 1000
tiv_ft  <- df$eTIV[df$Group == "FT"]  / 1000
tiv_vpt <- tiv_vpt[!is.na(tiv_vpt)]
tiv_ft  <- tiv_ft[!is.na(tiv_ft)]
tiv_t <- t.test(tiv_vpt, tiv_ft)
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "shared",
  Variable = "Total Intracranial Volume, cm\u00b3\u1d43",
  VPT = sprintf("%.0f \u00b1 %.0f", mean(tiv_vpt), sd(tiv_vpt)),
  FT  = sprintf("%.0f \u00b1 %.0f", mean(tiv_ft),  sd(tiv_ft)),
  p   = fmt_p(tiv_t$p.value)
)

# Relative motion
rm_row <- welch_row("Rel_Motion", "Relative Motion", "mm", digits = 2)
t1_rows[[length(t1_rows) + 1]] <- list(section = "shared",
                                       Variable = paste0(rm_row$label, "\u1d43"),
                                       VPT = rm_row$VPT, FT = rm_row$FT,
                                       p = rm_row$p)

# VPT clinical variables header
t1_rows[[length(t1_rows) + 1]] <- list(section = "vpt",
                                       Variable = "VPT clinical variables",
                                       VPT = "", FT = "", p = "")

# Birth weight (VPT only)
bw_vpt <- df$bw[df$Group == "VPT"]
bw_vpt <- bw_vpt[!is.na(bw_vpt)]
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "vpt_sub",
  Variable = "Birth weight, grams",
  VPT = sprintf("%.1f \u00b1 %.1f", mean(bw_vpt), sd(bw_vpt)),
  FT  = "\u2014",
  p   = "\u2014"
)

# Maternal social risk score
srs_vpt <- df$sriskscore[df$Group == "VPT"]
srs_vpt <- srs_vpt[!is.na(srs_vpt)]
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "vpt_sub",
  Variable = "Maternal social risk score",
  VPT = sprintf("%.1f \u00b1 %.1f", mean(srs_vpt), sd(srs_vpt)),
  FT  = "\u2014",
  p   = "\u2014"
)

# Global brain abnormality score (median [IQR])
gba_vpt <- df$globalbrainscore2[df$Group == "VPT"]
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "vpt_sub",
  Variable = "Global brain abnormality score\u1d9c",
  VPT = median_iqr(gba_vpt, digits = 1),
  FT  = "\u2014",
  p   = "\u2014"
)

# BPD grade
# Methods 2.1 says bpd2 is binary (1 = any BPD; 0 = none). Manuscript Table 1
# shows the underlying ordinal grade, carried in the bpdgrade column (0-3).
# If bpdgrade is absent we fall back to the binary bpd2 split.
if ("bpdgrade" %in% names(df)) {
  bpd_vpt <- df$bpdgrade[df$Group == "VPT"]
  t1_rows[[length(t1_rows) + 1]] <- list(
    section = "vpt", Variable = "BPD grade", VPT = "", FT = "", p = ""
  )
  for (g in c(0, 1, 2, 3)) {
    lbl <- if (g == 0) "  No BPD" else sprintf("  Grade %d", g)
    n_g <- sum(bpd_vpt == g, na.rm = TRUE)
    t1_rows[[length(t1_rows) + 1]] <- list(
      section = "vpt_sub", Variable = lbl,
      VPT = n_pct(n_g, n_vpt),
      FT  = "\u2014",
      p   = "\u2014"
    )
  }
} else {
  bpd_vpt <- df$bpd2[df$Group == "VPT"]
  t1_rows[[length(t1_rows) + 1]] <- list(
    section = "vpt", Variable = "BPD (any vs none)", VPT = "", FT = "", p = ""
  )
  t1_rows[[length(t1_rows) + 1]] <- list(
    section = "vpt_sub", Variable = "  No BPD",
    VPT = n_pct(sum(bpd_vpt == 0, na.rm = TRUE), n_vpt),
    FT = "\u2014", p = "\u2014"
  )
  t1_rows[[length(t1_rows) + 1]] <- list(
    section = "vpt_sub", Variable = "  Any BPD",
    VPT = n_pct(sum(bpd_vpt == 1, na.rm = TRUE), n_vpt),
    FT = "\u2014", p = "\u2014"
  )
}

# ROP
rop_vpt <- df$anyrop[df$Group == "VPT"]
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "vpt_sub",
  Variable = "ROP",
  VPT = n_pct(sum(rop_vpt == 1, na.rm = TRUE), n_vpt),
  FT  = "\u2014",
  p   = "\u2014"
)

# Sepsis
sep_vpt <- df$sepsis2[df$Group == "VPT"]
t1_rows[[length(t1_rows) + 1]] <- list(
  section = "vpt_sub",
  Variable = "Sepsis",
  VPT = n_pct(sum(sep_vpt == 1, na.rm = TRUE), n_vpt),
  FT  = "\u2014",
  p   = "\u2014"
)

# --- Assemble dataframe -----------------------------------------------------
table1_df <- data.frame(
  Variable  = sapply(t1_rows, `[[`, "Variable"),
  VPT       = sapply(t1_rows, `[[`, "VPT"),
  FT        = sapply(t1_rows, `[[`, "FT"),
  p_value   = sapply(t1_rows, `[[`, "p"),
  stringsAsFactors = FALSE
)
names(table1_df) <- c("",
                      sprintf("VPT (n = %d)", n_vpt),
                      sprintf("FT (n = %d)",  n_ft),
                      "p-value")

cat("Table 1 rows:\n")
print(table1_df, row.names = FALSE)

# --- Write CSV ---------------------------------------------------------------
out_path <- file.path(tables_out, "Table_1_participant_characteristics.csv")
write.csv(table1_df, out_path, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nWrote: %s\n", out_path))


# ============================================================================ 
# SECTION 4: TABLE 2 + FIGURE 2 - VPT vs FT GROUP DIFFERENCES
# ----------------------------------------------------------------------------
# Manuscript Table 2 layout (5 metric columns x 14 data rows):
#
#   Columns:  Density | Strength | GE | ACC | SW
#
#   Rows (per Unadjusted / Fully Adjusted block):
#     n (VPT/FT)
#     Test statistic               (t / F / chi^2 / beta depending on metric)
#     p-value
#     Variance explained
#     Functional variance explained
#     Effect direction             (VPT < FT  or  VPT > FT)
#     CI excludes zero             (e.g. "11-100%", "ns" - PFFR rows only)
#
# Density column sourced from results/density_group_difference/
# density_sensitivity_results.xlsx (Sheet "Sensitivity Analyses"):
#   Sample == "Full sample" + Covariates %in% c("Unadjusted (Welch's t)",
#                                              "+ Age + Sex + Motion + eTIV")
#
# PFFR metric columns sourced from 8 RData files in results/fda_group_differences/.
# Folder naming convention (TWO_fda_group_difference_runs.R):
#   {metric}_FDA_Group_11-100_fullsample                                (unadj)
#   {metric}_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample (full)
# where {metric} prefix is one of: strength, rand_norm_wei_GE,
# rand_norm_wei_ACC, rand_norm_wei_SW.
#
# Figure 2: 4-panel composite (A=Strength, B=GE, C=ACC, D=SW), each showing the
# Group beta(density) coefficient with 95% bootstrap CI from the fully-adjusted
# fit. Sourced from the same fully-adjusted RData files used for Table 2.
# ============================================================================ 

cat("\n========================================\n")
cat("Building Table 2 + Figure 2 (Group Differences)\n")
cat("========================================\n\n")

# --- Metric specs ----------------------------------------------------------
# Order = manuscript order. family_label drives test-statistic formatting:
#   gaussian -> "F"; scat (SW) -> "Chi^2".
# ve_suffix appends 'a' to VE percentages for scat metrics, matching the
# manuscript footnote convention.
group_specs <- list(
  Strength = list(
    folder_prefix  = "str",
    rdata_filename = "str_FDA_results.RData",
    panel_label    = "Strength",
    family_label   = "F",
    ve_suffix      = ""
  ),
  GE = list(
    folder_prefix  = "rand_norm_wei_GE",
    rdata_filename = "rand_norm_wei_GE_FDA_results.RData",
    panel_label    = "GE",
    family_label   = "F",
    ve_suffix      = ""
  ),
  ACC = list(
    folder_prefix  = "rand_norm_wei_ACC",
    rdata_filename = "rand_norm_wei_ACC_FDA_results.RData",
    panel_label    = "ACC",
    family_label   = "F",
    ve_suffix      = ""
  ),
  SW = list(
    folder_prefix  = "rand_norm_wei_SW",
    rdata_filename = "rand_norm_wei_SW_FDA_results.RData",
    panel_label    = "SW",
    family_label   = "Chi^2",
    ve_suffix      = "\u1d43"   # superscript 'a' footnote marker
  )
)

# Folder suffixes (predictor order MUST match TWO_fda_group_difference_runs.R)
suffix_unadj <- "Group_11-100_fullsample"
suffix_full  <- "Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample"

# --- Helpers (Section-4-scoped) ---------------------------------------------

# Format a PFFR test statistic with the family-appropriate symbol, e.g.
# "F = 11.95" or "chi^2 = 16.42^a".
fmt_pffr_stat <- function(stat, family_label, suffix = "") {
  symbol <- if (family_label == "Chi^2") "\u03c7\u00b2" else "F"
  paste0(symbol, " = ", sprintf("%.2f", stat), suffix)
}

# Format a variance-explained value (input on [0,1] or [0,100]; auto-detect).
# (Moved to Section 1 helpers so Section 5 can use it without sourcing Section 4.)

# Determine effect direction from the bootstrap Group beta(d) curve.
# Returns "VPT < FT" if beta is predominantly negative over densities where the
# 95% CI excludes zero, "VPT > FT" if positive. If the CI never excludes zero,
# falls back to the sign of the curve at its peak magnitude.
group_effect_direction <- function(beta_df) {
  if (nrow(beta_df) == 0) return(NA_character_)
  excludes <- (beta_df$ci_lower > 0) | (beta_df$ci_upper < 0)
  if (any(excludes, na.rm = TRUE)) {
    sig_betas <- beta_df$beta[excludes]
    sign_val <- sign(mean(sig_betas, na.rm = TRUE))
  } else {
    # No significant region: sign at peak |beta|
    peak_idx <- which.max(abs(beta_df$beta))
    sign_val <- sign(beta_df$beta[peak_idx])
  }
  if (is.na(sign_val) || sign_val == 0) return(NA_character_)
  if (sign_val < 0) "VPT < FT" else "VPT > FT"
}

# Compute Group n (VPT/FT) after metric-specific outlier exclusion.
# Reads predictor_df_clean (the post-exclusion design frame) from the RData env.
# Group is stored numerically (0=FT, 1=VPT) per Methods 2.1.
group_n_from_env <- function(env) {
  pdc <- env$predictor_df_clean
  if (is.null(pdc) || !"Group" %in% names(pdc)) {
    return(c(VPT = NA_integer_, FT = NA_integer_))
  }
  c(VPT = sum(pdc$Group == 1, na.rm = TRUE),
    FT  = sum(pdc$Group == 0, na.rm = TRUE))
}

# Extract one metric's "cells dictionary" for Table 2 (one config = unadjusted
# or fully adjusted). Also caches the fully-adjusted beta_df for Figure 2.
# Returns a list with keys: n, stat, p, ve, fve, direction, ci_range, beta_df.
extract_group_cells <- function(rdata_path, spec, is_fully_adjusted) {
  env <- load_rdata(rdata_path)
  
  sr <- env$summary_results
  if (is.null(sr)) stop(sprintf("No summary_results in %s", rdata_path))
  
  ns <- group_n_from_env(env)
  stat_val <- sr$Overall_Test_Statistic[1]
  p_val    <- sr$Overall_p[1]
  
  ve_val  <- sr$Variance_Explained[1]
  fve_val <- sr$Functional_VE[1]
  
  # beta(d) curve and CI exclusion range from bootstrap Group term
  beta_df <- tryCatch(
    extract_beta_curve(env, "Group", model_obj_name = "pffr_fit"),
    error = function(e) {
      warning(sprintf("Could not extract Group beta from %s: %s",
                      rdata_path, conditionMessage(e)))
      data.frame(density = numeric(0), beta = numeric(0),
                 ci_lower = numeric(0), ci_upper = numeric(0))
    }
  )
  
  ci_range <- ci_exclusion_range(beta_df)
  direction <- group_effect_direction(beta_df)
  
  list(
    n         = ns,
    stat      = stat_val,
    p         = p_val,
    ve        = ve_val,
    fve       = fve_val,
    direction = direction,
    ci_range  = ci_range,
    beta_df   = beta_df  # cached for Figure 2 if fully adjusted
  )
}

# --- Iterate: build PFFR cells + collect Figure 2 envs ----------------------

pffr_cells <- list()    # pffr_cells[[mkey]] = list(unadj=..., full=...)
fig2_betas <- list()    # fig2_betas[[mkey]] = beta_df for fully-adjusted Group

for (mkey in names(group_specs)) {
  spec <- group_specs[[mkey]]
  cat(sprintf("[%s]\n", mkey))
  
  unadj_path <- file.path(fda_group_root,
                          sprintf("%s_FDA_%s", spec$folder_prefix, suffix_unadj),
                          spec$rdata_filename)
  full_path  <- file.path(fda_group_root,
                          sprintf("%s_FDA_%s", spec$folder_prefix, suffix_full),
                          spec$rdata_filename)
  
  if (require_input(unadj_path, sprintf("%s unadjusted RData", mkey))) {
    unadj_cells <- extract_group_cells(unadj_path, spec, is_fully_adjusted = FALSE)
    cat(sprintf("  Unadj  stat=%.2f p=%.4g  n(VPT/FT)=%d/%d  CI=%s\n",
                unadj_cells$stat, unadj_cells$p,
                unadj_cells$n["VPT"], unadj_cells$n["FT"],
                unadj_cells$ci_range))
  } else {
    unadj_cells <- NULL
  }
  
  if (require_input(full_path, sprintf("%s fully-adjusted RData", mkey))) {
    full_cells <- extract_group_cells(full_path, spec, is_fully_adjusted = TRUE)
    cat(sprintf("  Adj    stat=%.2f p=%.4g  n(VPT/FT)=%d/%d  CI=%s\n",
                full_cells$stat, full_cells$p,
                full_cells$n["VPT"], full_cells$n["FT"],
                full_cells$ci_range))
    fig2_betas[[mkey]] <- full_cells$beta_df
  } else {
    full_cells <- NULL
  }
  
  pffr_cells[[mkey]] <- list(unadj = unadj_cells, full = full_cells)
}

# --- Density column from density_sensitivity_results.xlsx -------------------
# Sheet "Sensitivity Analyses". We pull rows where Sample == "Full sample" and
# Covariates %in% c("Unadjusted (Welch's t)", "+ Age + Sex + Motion + eTIV").
# Unadjusted row reports |t| from Welch's; fully adjusted reports beta.

density_cells <- list(unadj = NULL, full = NULL)
if (require_input(density_xlsx, "density results xlsx")) {
  dens_res <- as.data.frame(read_excel(density_xlsx, sheet = "Sensitivity Analyses"))
  
  # Robust column lookup (some R versions strip parentheses / replace spaces)
  pick_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) {
      stop(sprintf("Density xlsx missing any of: %s",
                   paste(candidates, collapse = ", ")))
    }
    hit[1]
  }
  col_sample  <- pick_col(dens_res, c("Sample"))
  col_covs    <- pick_col(dens_res, c("Covariates"))
  col_beta    <- pick_col(dens_res, c("Beta (Group)", "Beta"))
  col_t       <- pick_col(dens_res, c("t"))
  col_p       <- pick_col(dens_res, c("p"))
  col_n_vpt   <- pick_col(dens_res, c("N (VPT)", "N_VPT"))
  col_n_ft    <- pick_col(dens_res, c("N (FT)",  "N_FT"))
  
  full_sample <- dens_res[dens_res[[col_sample]] == "Full sample", , drop = FALSE]
  
  unadj_row <- full_sample[full_sample[[col_covs]] == "Unadjusted (Welch's t)",
                           , drop = FALSE]
  adj_row   <- full_sample[full_sample[[col_covs]] == "+ Age + Sex + Motion + eTIV",
                           , drop = FALSE]
  
  if (nrow(unadj_row) == 1) {
    t_val <- abs(unadj_row[[col_t]][1])  # report magnitude only, sign in direction row
    density_cells$unadj <- list(
      n         = c(VPT = unadj_row[[col_n_vpt]][1], FT = unadj_row[[col_n_ft]][1]),
      stat_str  = sprintf("t = %.2f", t_val),
      p         = unadj_row[[col_p]][1],
      ve        = NA_real_,
      fve       = NA_real_,
      direction = if (unadj_row[[col_t]][1] < 0) "VPT > FT" else "VPT < FT"
    )
    # Use fmt_p_relop for the diagnostic so character p ("<0.001") prints fine.
    cat(sprintf("Density unadj:  t=%.2f, p=%s  n(VPT/FT)=%d/%d\n",
                t_val, fmt_p_relop(density_cells$unadj$p),
                density_cells$unadj$n["VPT"], density_cells$unadj$n["FT"]))
  }
  if (nrow(adj_row) == 1) {
    beta_val <- adj_row[[col_beta]][1]
    density_cells$full <- list(
      n         = c(VPT = adj_row[[col_n_vpt]][1], FT = adj_row[[col_n_ft]][1]),
      stat_str  = sprintf("\u03b2 = %s%.3f", if (beta_val < 0) "\u2212" else "",
                          abs(beta_val)),
      p         = adj_row[[col_p]][1],
      ve        = NA_real_,
      fve       = NA_real_,
      direction = if (beta_val < 0) "VPT < FT" else "VPT > FT"
    )
    cat(sprintf("Density adj:    beta=%.3f, p=%s\n",
                beta_val, fmt_p_relop(density_cells$full$p)))
  }
}

# --- Assemble Table 2 dataframe ---------------------------------------------
# Each "cell builder" returns the formatted string for one row x metric.

cell_n        <- function(c) if (is.null(c)) "PENDING" else sprintf("%d/%d", c$n["VPT"], c$n["FT"])
cell_stat_pffr<- function(c, spec) if (is.null(c)) "PENDING" else fmt_pffr_stat(c$stat, spec$family_label, spec$ve_suffix)
cell_p        <- function(c) if (is.null(c)) "PENDING" else fmt_p_relop(c$p)
cell_ve       <- function(c, spec) if (is.null(c)) "PENDING" else fmt_ve(c$ve, spec$ve_suffix)
cell_fve      <- function(c, spec) if (is.null(c)) "PENDING" else fmt_ve(c$fve, spec$ve_suffix)
cell_dir      <- function(c) if (is.null(c)) "PENDING" else c$direction
cell_ci       <- function(c) if (is.null(c)) "PENDING" else c$ci_range

# Density-only formatters (no PFFR statistic, no VE/FVE/CI)
cell_stat_den <- function(c) if (is.null(c)) "PENDING" else c$stat_str

row_dash <- function() "\u2014"

mkeys <- names(group_specs)
table2_rows <- list()

# Header row: section "Unadjusted"
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Unadjusted", Density = "",
  setNames(rep("", length(mkeys)), mkeys)
)

# Unadjusted block
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "n (VPT/FT)",
  Density  = cell_n(density_cells$unadj),
  sapply(mkeys, function(m) cell_n(pffr_cells[[m]]$unadj))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Test statistic",
  Density  = cell_stat_den(density_cells$unadj),
  sapply(mkeys, function(m) cell_stat_pffr(pffr_cells[[m]]$unadj, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "p-value",
  Density  = cell_p(density_cells$unadj),
  sapply(mkeys, function(m) cell_p(pffr_cells[[m]]$unadj))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Variance explained",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_ve(pffr_cells[[m]]$unadj, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Functional variance explained",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_fve(pffr_cells[[m]]$unadj, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Effect direction",
  Density  = if (is.null(density_cells$unadj)) "PENDING" else density_cells$unadj$direction,
  sapply(mkeys, function(m) cell_dir(pffr_cells[[m]]$unadj))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "CI excludes zero",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_ci(pffr_cells[[m]]$unadj))
)

# Header row: section "Fully Adjusted"
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Fully Adjusted", Density = "",
  setNames(rep("", length(mkeys)), mkeys)
)

# Fully-adjusted block
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "n (VPT/FT)",
  Density  = cell_n(density_cells$full),
  sapply(mkeys, function(m) cell_n(pffr_cells[[m]]$full))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Test statistic",
  Density  = cell_stat_den(density_cells$full),
  sapply(mkeys, function(m) cell_stat_pffr(pffr_cells[[m]]$full, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "p-value",
  Density  = cell_p(density_cells$full),
  sapply(mkeys, function(m) cell_p(pffr_cells[[m]]$full))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Variance explained",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_ve(pffr_cells[[m]]$full, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Functional variance explained",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_fve(pffr_cells[[m]]$full, group_specs[[m]]))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "Effect direction",
  Density  = if (is.null(density_cells$full)) "PENDING" else density_cells$full$direction,
  sapply(mkeys, function(m) cell_dir(pffr_cells[[m]]$full))
)
table2_rows[[length(table2_rows) + 1]] <- c(
  Variable = "CI excludes zero",
  Density  = row_dash(),
  sapply(mkeys, function(m) cell_ci(pffr_cells[[m]]$full))
)

table2_df <- as.data.frame(do.call(rbind, table2_rows), stringsAsFactors = FALSE)
names(table2_df) <- c("", "Density", "Strength", "GE", "ACC", "SW")

cat("\nTable 2:\n")
print(table2_df, row.names = FALSE)

out_path <- file.path(tables_out, "Table_2_group_differences.csv")
write.csv(table2_df, out_path, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nWrote: %s\n", out_path))


# --- Figure 2: 4-panel composite of fully-adjusted Group beta(d) -----------

build_fig2_panel <- function(beta_df, panel_letter, panel_title) {
  if (is.null(beta_df) || nrow(beta_df) == 0) {
    # Placeholder panel when RData not yet available
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = sprintf("%s. %s\n(PENDING)", panel_letter, panel_title)) +
             theme_void())
  }
  
  # beta(density) ribbon panel: point estimate with bootstrap 95% CI band.
  ggplot(beta_df, aes(x = density, y = beta)) +
    geom_hline(yintercept = 0, lty = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.3, fill = "steelblue") +
    geom_line(linewidth = 1, color = "darkblue") +
    ylab(expression(beta(d))) +
    xlab("Density (%)") +
    ggtitle(sprintf("%s. %s", panel_letter, panel_title)) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}

panel_letters <- c(Strength = "A", GE = "B", ACC = "C", SW = "D")
fig2_panels <- list()
for (m in mkeys) {
  fig2_panels[[m]] <- build_fig2_panel(
    fig2_betas[[m]], panel_letters[m], group_specs[[m]]$panel_label
  )
}
fig2 <- (fig2_panels$Strength | fig2_panels$GE) /
  (fig2_panels$ACC      | fig2_panels$SW)

save_fig(fig2, "Figure_2_group_difference_betas",
         width_in = 10, height_in = 7)


# ============================================================================ 
# SECTION 5: TABLE 3 + FIGURES 3-5 - STABILITY SELECTION
# ----------------------------------------------------------------------------
# Manuscript Table 3 layout (4 metric columns, ordering = stably-selected first):
#   Columns:  ACC | SW | GE | Strength
#   Row sections:
#     Stability selection
#       Top exposure (frequency)        e.g. "GBA (94%)"
#       2nd exposure (frequency)        e.g. "BPD (58%)"
#     Full-sample model
#       Test statistic                  F = ... or chi^2 = ...^a
#       p
#       Variance explained
#       Functional variance explained
#     Selected exposure contribution
#       Test statistic (vs forced only)
#       p (vs forced only)
#       Functional variance explained increment
#       Beta direction                  Positive / Negative
#       CI excludes zero                density range "11-100%" / "15.5-68.5%"
# For metrics without a stably-selected exposure (GE, Strength), the model
# and contribution rows show dashes; only the two exposure rows are filled.
#
# Inputs: 4 RData files at
#   results/fda_stability_selection/{metric}_stabsel_11-100/
#       {metric}_stabsel_results.RData
# where {metric} is: str, rand_norm_wei_GE, rand_norm_wei_ACC, rand_norm_wei_SW.
#
# Figure layouts:
#   Figure 3: ACC — 2 panels (A: stability frequencies, B: GBA beta(d))
#   Figure 4: SW  — 2 panels (A: stability frequencies, B: ROP beta(d))
#   Figure 5: GE | Strength — 2 panels (both stability frequencies only)
# ============================================================================ 

cat("\n========================================\n")
cat("Building Table 3 + Figures 3, 4, 5 (Stability Selection)\n")
cat("========================================\n\n")

# --- Metric specs (Table 3 column order: ACC, SW, GE, Strength) -------------
# folder_subdir is relative to fda_stabsel_root; rdata_filename comes from the
# pipeline's save() call (see FOUR_fda_stability_selection_pipeline.R line 1009).
stabsel_specs <- list(
  ACC = list(
    folder_subdir  = "rand_norm_wei_ACC_stabsel_11-100",
    rdata_filename = "rand_norm_wei_ACC_stabsel_results.RData",
    family_label   = "F",
    ve_suffix      = "",
    metric_label   = "ACC",
    has_figure     = TRUE,
    figure_number  = 3,
    figure_basename = "Figure_3_ACC_stability_selection"
  ),
  SW = list(
    folder_subdir  = "rand_norm_wei_SW_stabsel_11-100",
    rdata_filename = "rand_norm_wei_SW_stabsel_results.RData",
    family_label   = "Chi^2",
    ve_suffix      = "\u1d43",   # superscript 'a' (scat footnote marker)
    metric_label   = "SW",
    has_figure     = TRUE,
    figure_number  = 4,
    figure_basename = "Figure_4_SW_stability_selection"
  ),
  GE = list(
    folder_subdir  = "rand_norm_wei_GE_stabsel_11-100",
    rdata_filename = "rand_norm_wei_GE_stabsel_results.RData",
    family_label   = "F",
    ve_suffix      = "",
    metric_label   = "GE",
    has_figure     = FALSE,   # combined into Figure 5 panel A
    figure_number  = NA,
    figure_basename = NA
  ),
  Strength = list(
    folder_subdir  = "str_stabsel_11-100",
    rdata_filename = "str_stabsel_results.RData",
    family_label   = "F",
    ve_suffix      = "",
    metric_label   = "Strength",
    has_figure     = FALSE,   # combined into Figure 5 panel B
    figure_number  = NA,
    figure_basename = NA
  )
)

# --- Exposure display labels (manuscript convention) ------------------------
exposure_labels <- c(
  bpd2              = "BPD",
  bw_z              = "BWZ",
  ga                = "GA",
  globalbrainscore2 = "GBA",
  anyrop            = "ROP",
  sepsis2           = "Sepsis"
)

# Stability threshold (Methods 2.7.4 - matches FOUR_fda_stability_selection_pipeline.R)
stability_threshold <- 0.70

# --- Helpers (Section-5-scoped) ---------------------------------------------

# Convert "globalbrainscore2" or other raw exposure name to display label.
label_exposure <- function(raw_name) {
  if (raw_name %in% names(exposure_labels)) exposure_labels[[raw_name]]
  else raw_name
}

# Format frequency as integer percent ("94%").
fmt_pct_int <- function(x) {
  if (is.na(x)) return("\u2014")
  sprintf("%d%%", as.integer(round(x * 100)))
}

# Pull the top-n exposures by selection frequency. Returns a data.frame with
# columns: exposure_raw, exposure_display, frequency.
top_exposures <- function(selection_freq, n = 2) {
  ord <- order(selection_freq, decreasing = TRUE)
  raw  <- names(selection_freq)[ord][seq_len(n)]
  freq <- selection_freq[ord][seq_len(n)]
  data.frame(
    exposure_raw     = raw,
    exposure_display = vapply(raw, label_exposure, ""),
    frequency        = as.numeric(freq),
    stringsAsFactors = FALSE
  )
}

# Format "Exposure (xx%)" for Table 3 cells.
fmt_exposure_freq <- function(display, freq) {
  sprintf("%s (%s)", display, fmt_pct_int(freq))
}

# Determine effect direction from bootstrap CI of the stably-selected exposure.
# Uses the same logic as group_effect_direction: sign at densities where CI
# excludes zero; if none, sign at peak |beta|.
exposure_direction <- function(beta_df) {
  if (nrow(beta_df) == 0) return(NA_character_)
  excludes <- (beta_df$ci_lower > 0) | (beta_df$ci_upper < 0)
  if (any(excludes, na.rm = TRUE)) {
    sig_betas <- beta_df$beta[excludes]
    sign_val <- sign(mean(sig_betas, na.rm = TRUE))
  } else {
    peak_idx <- which.max(abs(beta_df$beta))
    sign_val <- sign(beta_df$beta[peak_idx])
  }
  if (is.na(sign_val) || sign_val == 0) return(NA_character_)
  if (sign_val < 0) "Negative" else "Positive"
}

# Format Table 3 test-statistic cell. For SW (scat), uses chi^2 + footnote 'a';
# for gaussian, uses F.
fmt_stabsel_stat <- function(stat, spec) {
  if (is.na(stat)) return("\u2014")
  symbol <- if (spec$family_label == "Chi^2") "\u03c7\u00b2" else "F"
  paste0(symbol, " = ", sprintf("%.2f", stat), spec$ve_suffix)
}

# Selection-frequency vertical bar plot for one metric.
build_stabsel_freq_panel <- function(selection_freq, panel_letter, panel_title) {
  if (is.null(selection_freq) || length(selection_freq) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = sprintf("%s. %s\n(PENDING)", panel_letter, panel_title)) +
             theme_void())
  }
  
  freq_df <- data.frame(
    exposure_raw     = names(selection_freq),
    exposure_display = vapply(names(selection_freq), label_exposure, ""),
    frequency_pct    = as.numeric(selection_freq) * 100,
    stringsAsFactors = FALSE
  )
  freq_df <- freq_df[order(freq_df$frequency_pct, decreasing = TRUE), ]
  freq_df$exposure_display <- factor(freq_df$exposure_display,
                                     levels = freq_df$exposure_display)
  freq_df$above <- freq_df$frequency_pct >= stability_threshold * 100
  freq_df$label <- sprintf("%d%%", as.integer(round(freq_df$frequency_pct)))
  
  ggplot(freq_df, aes(x = exposure_display, y = frequency_pct,
                      fill = above)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = stability_threshold * 100,
               lty = "dashed", color = "grey50", linewidth = 0.6) +
    geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("FALSE" = "#5b89c7", "TRUE" = "#1a237e"),
                      guide = "none") +
    scale_y_continuous(limits = c(0, 105),
                       breaks = c(0, 25, 50, 75, 100),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Selection Frequency (%)",
         title = sprintf("%s. %s", panel_letter, panel_title)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 13),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 11),
          axis.ticks.x = element_blank())
}

# beta(d) panel for a stably selected exposure.
build_stabsel_beta_panel <- function(beta_df, panel_letter, panel_title) {
  if (is.null(beta_df) || nrow(beta_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = sprintf("%s. %s\n(PENDING)", panel_letter, panel_title)) +
             theme_void())
  }
  ggplot(beta_df, aes(x = density, y = beta)) +
    geom_hline(yintercept = 0, lty = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.3, fill = "steelblue") +
    geom_line(linewidth = 1, color = "darkblue") +
    ylab(expression(beta(d))) +
    xlab("Density (%)") +
    ggtitle(sprintf("%s. %s", panel_letter, panel_title)) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}

# --- Iterate over metrics: extract cells, cache panels ---------------------
# stabsel_cells[[mkey]] = list of named character cells + diagnostic fields
# fig_panels[[mkey]]   = list(freq_panel = ggplot, beta_panel = ggplot or NULL)

stabsel_cells <- list()
fig_panels    <- list()

mkeys_stabsel <- names(stabsel_specs)

for (mkey in mkeys_stabsel) {
  spec <- stabsel_specs[[mkey]]
  cat(sprintf("[%s]\n", mkey))
  
  rdata_path <- file.path(fda_stabsel_root, spec$folder_subdir, spec$rdata_filename)
  
  if (!require_input(rdata_path, sprintf("%s stabsel RData", mkey))) {
    # Dev mode: emit pending placeholder cells, skip figure generation
    stabsel_cells[[mkey]] <- NULL
    fig_panels[[mkey]]    <- list(freq_panel = NULL, beta_panel = NULL)
    next
  }
  
  env <- load_rdata(rdata_path)
  
  # Selection frequency: pipeline stores as a named numeric vector. Confirm.
  selfreq <- env$selection_freq
  if (is.null(selfreq)) {
    if (STRICT) stop(sprintf("No selection_freq in %s", rdata_path))
    cat("  [skip] no selection_freq present\n")
    stabsel_cells[[mkey]] <- NULL
    next
  }
  
  # Top two exposures (for Table 3 rows + figure context)
  top2 <- top_exposures(selfreq, n = 2)
  top1_cell <- fmt_exposure_freq(top2$exposure_display[1], top2$frequency[1])
  top2_cell <- fmt_exposure_freq(top2$exposure_display[2], top2$frequency[2])
  cat(sprintf("  Top exposure: %s\n", top1_cell))
  cat(sprintf("  2nd exposure: %s\n", top2_cell))
  
  # Stable exposures (manuscript convention: at most one per metric in main results)
  stable <- env$stable_exposures
  if (is.null(stable)) stable <- character(0)
  has_stable <- length(stable) > 0
  
  # Cells that only exist if at least one exposure stably selected
  if (has_stable) {
    # Overall full-sample model stats
    overall_stat <- env$overall_test_stat
    overall_p_val <- env$overall_p
    ve_val  <- env$ve_full
    fve_val <- env$functional_ve
    
    # Selected exposure contribution stats
    lr_stat <- env$lr_stat_vs_forced
    lr_p    <- env$lr_p_vs_forced
    fve_inc <- env$functional_ve_increment
    
    # Beta(d) curve for the stably-selected exposure. 
    sel_exp <- stable[1]   # In the main analysis only one stable per metric
    if (is.null(env$bootstrap_coefs)) {
      cat(sprintf("  [warn] %s has stable exposure %s but bootstrap_coefs is NULL\n",
                  mkey, sel_exp))
      cat("         direction + CI exclusion range will be em-dashes.\n")
      beta_df   <- NULL
      direction <- NA_character_
      ci_range  <- NA_character_
    } else {
      beta_df <- tryCatch(
        extract_beta_curve(env, sel_exp, model_obj_name = "pffr_final"),
        error = function(e) {
          warning(sprintf("Could not extract %s beta from %s: %s",
                          sel_exp, rdata_path, conditionMessage(e)))
          data.frame(density = numeric(0), beta = numeric(0),
                     ci_lower = numeric(0), ci_upper = numeric(0))
        }
      )
      direction <- exposure_direction(beta_df)
      ci_range  <- ci_exclusion_range(beta_df)
      
      cat(sprintf("  Stably selected: %s -> direction %s, CI %s\n",
                  label_exposure(sel_exp), direction, ci_range))
    }
  } else {
    overall_stat <- NA_real_; overall_p_val <- NA_real_
    ve_val <- NA_real_; fve_val <- NA_real_
    lr_stat <- NA_real_; lr_p <- NA_real_; fve_inc <- NA_real_
    sel_exp <- NA_character_; beta_df <- NULL
    direction <- NA_character_; ci_range <- NA_character_
    cat("  No stably-selected exposure.\n")
  }
  
  # Assemble formatted cells
  cells <- list(
    top_exposure         = top1_cell,
    second_exposure      = top2_cell,
    test_statistic       = if (has_stable) fmt_stabsel_stat(overall_stat, spec) else "\u2014",
    p_overall            = if (has_stable) fmt_p_relop(overall_p_val) else "\u2014",
    variance_explained   = if (has_stable) fmt_ve(ve_val, spec$ve_suffix) else "\u2014",
    functional_ve        = if (has_stable) fmt_ve(fve_val, spec$ve_suffix) else "\u2014",
    test_statistic_vs_forced = if (has_stable) fmt_stabsel_stat(lr_stat, spec) else "\u2014",
    p_vs_forced          = if (has_stable) fmt_p_relop(lr_p) else "\u2014",
    fve_increment        = if (has_stable) fmt_ve(fve_inc, spec$ve_suffix) else "\u2014",
    beta_direction       = if (has_stable) direction else "\u2014",
    ci_excludes_zero     = if (has_stable) ci_range else "\u2014"
  )
  stabsel_cells[[mkey]] <- cells
  
  # Build figure panels
  # Title format matches manuscript: "ACC stability selection" (Fig 3A,4A,5)
  freq_title <- sprintf("%s stability selection", spec$metric_label)
  freq_panel <- build_stabsel_freq_panel(selfreq, "A", freq_title)
  if (has_stable) {
    beta_title <- sprintf("%s effect on %s",
                          label_exposure(sel_exp), spec$metric_label)
    beta_panel <- build_stabsel_beta_panel(beta_df, "B", beta_title)
  } else {
    beta_panel <- NULL
  }
  fig_panels[[mkey]] <- list(freq_panel = freq_panel, beta_panel = beta_panel)
  
  cat("\n")
}

# --- Assemble Table 3 dataframe ---------------------------------------------
# Manuscript ordering: ACC, SW, GE, Strength.

cell_or_pending <- function(cells, key) {
  if (is.null(cells)) return("PENDING")
  cells[[key]]
}

mkey_order_t3 <- c("ACC", "SW", "GE", "Strength")
metric_n <- list(    
  ACC      = 171,
  SW       = 169,
  GE       = 170,
  Strength = 169
)

table3_rows <- list()

# Header for "Stability selection" section (italic block label)
table3_rows[[length(table3_rows) + 1]] <- c(
  Variable = "Stability selection",
  setNames(rep("", length(mkey_order_t3)), mkey_order_t3)
)
table3_rows[[length(table3_rows) + 1]] <- c(
  Variable = "Top exposure (frequency)",
  setNames(vapply(mkey_order_t3,
                  function(m) cell_or_pending(stabsel_cells[[m]], "top_exposure"),
                  ""),
           mkey_order_t3)
)
table3_rows[[length(table3_rows) + 1]] <- c(
  Variable = "2nd exposure (frequency)",
  setNames(vapply(mkey_order_t3,
                  function(m) cell_or_pending(stabsel_cells[[m]], "second_exposure"),
                  ""),
           mkey_order_t3)
)

# Header for "Full-sample model"
table3_rows[[length(table3_rows) + 1]] <- c(
  Variable = "Full-sample model",
  setNames(rep("", length(mkey_order_t3)), mkey_order_t3)
)
for (entry in list(
  list(label = "Test statistic",                key = "test_statistic"),
  list(label = "p",                              key = "p_overall"),
  list(label = "Variance explained",             key = "variance_explained"),
  list(label = "Functional variance explained",  key = "functional_ve"))) {
  table3_rows[[length(table3_rows) + 1]] <- c(
    Variable = entry$label,
    setNames(vapply(mkey_order_t3,
                    function(m) cell_or_pending(stabsel_cells[[m]], entry$key),
                    ""),
             mkey_order_t3)
  )
}

# Header for "Selected exposure contribution"
table3_rows[[length(table3_rows) + 1]] <- c(
  Variable = "Selected exposure contribution",
  setNames(rep("", length(mkey_order_t3)), mkey_order_t3)
)
for (entry in list(
  list(label = "Test statistic (vs forced covariate model only)",
       key = "test_statistic_vs_forced"),
  list(label = "p (vs forced covariate model only)",
       key = "p_vs_forced"),
  list(label = "Functional variance explained increment",
       key = "fve_increment"),
  list(label = "Beta direction",                key = "beta_direction"),
  list(label = "CI excludes zero",              key = "ci_excludes_zero"))) {
  table3_rows[[length(table3_rows) + 1]] <- c(
    Variable = entry$label,
    setNames(vapply(mkey_order_t3,
                    function(m) cell_or_pending(stabsel_cells[[m]], entry$key),
                    ""),
             mkey_order_t3)
  )
}

table3_df <- as.data.frame(do.call(rbind, table3_rows), stringsAsFactors = FALSE)
names(table3_df) <- c("",
                      sprintf("ACC (n=%d)", metric_n$ACC),
                      sprintf("SW (n=%d)",  metric_n$SW),
                      sprintf("GE (n=%d)",  metric_n$GE),
                      sprintf("Strength (n=%d)", metric_n$Strength))

cat("\nTable 3:\n")
print(table3_df, row.names = FALSE)

out_path <- file.path(tables_out, "Table_3_stability_selection.csv")
write.csv(table3_df, out_path, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nWrote: %s\n", out_path))

# --- Build Figure 3 (ACC stab freq + GBA beta(d)) ---------------------------

if (!is.null(fig_panels$ACC$freq_panel)) {
  if (!is.null(fig_panels$ACC$beta_panel)) {
    fig3 <- fig_panels$ACC$freq_panel | fig_panels$ACC$beta_panel
  } else {
    fig3 <- fig_panels$ACC$freq_panel
  }
  save_fig(fig3, stabsel_specs$ACC$figure_basename,
           width_in = 10, height_in = 4.5)
}

# --- Build Figure 4 (SW stab freq + ROP beta(d)) ----------------------------

if (!is.null(fig_panels$SW$freq_panel)) {
  if (!is.null(fig_panels$SW$beta_panel)) {
    fig4 <- fig_panels$SW$freq_panel | fig_panels$SW$beta_panel
  } else {
    fig4 <- fig_panels$SW$freq_panel
  }
  save_fig(fig4, stabsel_specs$SW$figure_basename,
           width_in = 10, height_in = 4.5)
}

# --- Build Figure 5 (GE stab freq | Strength stab freq) ---------------------

if (!is.null(fig_panels$GE$freq_panel) && !is.null(fig_panels$Strength$freq_panel)) {
  ge_panel <- fig_panels$GE$freq_panel +
    ggtitle("A. GE stability selection")
  str_panel <- fig_panels$Strength$freq_panel +
    ggtitle("B. Strength stability selection")
  fig5 <- ge_panel | str_panel
  save_fig(fig5, "Figure_5_GE_strength_stability",
           width_in = 10, height_in = 4.5)
} else if (!is.null(fig_panels$GE$freq_panel)) {
  save_fig(fig_panels$GE$freq_panel + ggtitle("A. GE stability selection"),
           "Figure_5_GE_strength_stability", width_in = 10, height_in = 4.5)
} else if (!is.null(fig_panels$Strength$freq_panel)) {
  save_fig(fig_panels$Strength$freq_panel + ggtitle("B. Strength stability selection"),
           "Figure_5_GE_strength_stability", width_in = 10, height_in = 4.5)
}