# ============================================================
# ONE_density_group_difference.R
#
# PURPOSE: Test group differences in unthresholded connectome density
#          between VPT and FT (Methods 2.7.1).
#
#          Primary test:
#            - Welch's t-test on raw density (VPT vs FT)
#
#          Sensitivity analyses:
#            - Five ANCOVAs, covarying individually and simultaneously
#              for corrected age at MRI, sex, relative motion, and TIV.
#              Partial eta-squared from Type II sums of squares.
#            - Low-density outlier robustness: identify VPT participants
#              with density > 2 SD below the VPT mean; re-run both
#              Welch's t and ANCOVAs on the trimmed cohort.
#
#          Produces:
#            1. Excel workbook of statistical results (3 sheets:
#               Results, Descriptives, Outliers)
#            2. Figure 1 (main paper) -- boxplot + forest plot
#
# INPUTS:
#   - IN_FILE: data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx
#
# OUTPUTS:
#   - OUT_XLSX: results/density_group_difference/density_sensitivity_results.xlsx
#   - OUT_FIG:  figures/main/Figure_1_density_group_comparison.png
#
# REQUIRES: R packages readxl, car, openxlsx, ggplot2, cowplot
# ============================================================

# ============================================================
# CONFIG -- edit before running
# ============================================================
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

IN_FILE  <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
OUT_XLSX <- file.path(repo_root, "results/density_group_difference/density_sensitivity_results.xlsx")
OUT_FIG  <- file.path(repo_root, "figures/main/Figure_1_density_group_comparison.png")
FORCE    <- FALSE   # set TRUE to overwrite existing outputs

# ============================================================
# Logic below -- do not edit
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(car)
  library(openxlsx)
  library(ggplot2)
  library(cowplot)
})

# --- Idempotency ---
if (file.exists(OUT_XLSX) && file.exists(OUT_FIG) && !FORCE) {
  cat("Outputs already exist:\n  ", OUT_XLSX, "\n  ", OUT_FIG, "\n",
      "Set FORCE=TRUE to overwrite. Exiting.\n", sep = "")
  quit(save = "no")
}

# --- Ensure output dirs exist ---
for (p in c(OUT_XLSX, OUT_FIG)) {
  d <- dirname(p)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ============================================================
# Load & prep
# ============================================================
df <- read_excel(IN_FILE)
df$Group   <- factor(df$Group, levels = c(0, 1), labels = c("FT", "VPT"))
df$density <- df$den_100.00
df$age5y   <- df$age_at_5y_mri
df$FD      <- df$Rel_Motion

# ============================================================
# 2 SD outlier identification (within VPT)
# ============================================================
vpt_dens   <- df$density[df$Group == "VPT"]
cutoff_2sd <- mean(vpt_dens) - 2 * sd(vpt_dens)
df$outlier  <- df$Group == "VPT" & df$density < cutoff_2sd
outlier_ids <- df$ID[df$outlier]
df_clean    <- df[!df$outlier, ]

cat("2 SD cutoff:", round(cutoff_2sd, 4),
    "| N outliers:", sum(df$outlier),
    "| IDs:", paste(outlier_ids, collapse = ", "), "\n\n")

# ============================================================
# Analysis function
# ============================================================
run_sensitivity <- function(data, label) {
  n_vpt <- sum(data$Group == "VPT")
  n_FT  <- sum(data$Group == "FT")
  
  # Welch's t-test (unadjusted)
  vpt_d <- data$density[data$Group == "VPT"]
  FT_d  <- data$density[data$Group == "FT"]
  tt    <- t.test(density ~ Group, data = data)
  pool  <- sqrt(((sd(vpt_d)^2 * (length(vpt_d)-1)) +
                   (sd(FT_d)^2  * (length(FT_d)-1))) /
                  (length(vpt_d) + length(FT_d) - 2))
  d     <- (mean(vpt_d) - mean(FT_d)) / pool
  
  unadj <- data.frame(
    Sample = label, N_VPT = n_vpt, N_FT = n_FT,
    Covariates = "Unadjusted (Welch's t)",
    Beta = NA_real_, t = round(tt$statistic[[1]], 2),
    p = tt$p.value, Partial_eta2 = NA_real_,
    Cohens_d = round(d, 3), R2 = NA_real_
  )
  
  # ANCOVA models: each covariate individually, then fully adjusted
  specs <- list(
    "+ Age at MRI"                  = density ~ Group + age5y,
    "+ Sex"                         = density ~ Group + sex,
    "+ Relative motion"             = density ~ Group + FD,
    "+ eTIV"                        = density ~ Group + eTIV,
    "+ Age + Sex + Motion + eTIV"   = density ~ Group + age5y + sex + FD + eTIV
  )
  
  ancova <- do.call(rbind, lapply(names(specs), function(nm) {
    m    <- lm(specs[[nm]], data = data)
    s    <- summary(m)
    aov2 <- Anova(m, type = 2)
    
    idx     <- grep("Group", rownames(s$coefficients))
    beta    <- s$coefficients[idx, "Estimate"]
    tval    <- s$coefficients[idx, "t value"]
    pval    <- s$coefficients[idx, "Pr(>|t|)"]
    ss_grp  <- aov2["Group", "Sum Sq"]
    ss_res  <- aov2["Residuals", "Sum Sq"]
    peta2   <- ss_grp / (ss_grp + ss_res)
    
    data.frame(
      Sample = label, N_VPT = n_vpt, N_FT = n_FT,
      Covariates = nm, Beta = round(beta, 4),
      t = round(tval, 2), p = pval,
      Partial_eta2 = round(peta2, 3),
      Cohens_d = NA_real_, R2 = round(s$r.squared, 3)
    )
  }))
  
  rbind(unadj, ancova)
}

# ============================================================
# Run sensitivity analyses
# ============================================================
tab <- rbind(
  run_sensitivity(df,       "Full sample"),
  run_sensitivity(df_clean, "Excl. outliers (2 SD)")
)
tab$p_formatted <- ifelse(tab$p < 0.001, "<.001", sprintf("%.4f", tab$p))

cat("Sensitivity results:\n")
print(tab[, c("Sample","Covariates","N_VPT","N_FT",
              "Beta","t","p_formatted","Partial_eta2","Cohens_d","R2")],
      row.names = FALSE)
cat("\n")

# ============================================================
# Write Excel workbook
# ============================================================
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "Bottom")
cs <- createStyle(halign = "center")

# Sheet 1: Results
addWorksheet(wb, "Sensitivity Analyses")
out_tab <- tab[, c("Sample","Covariates","N_VPT","N_FT",
                   "Beta","t","p_formatted","Partial_eta2","Cohens_d","R2")]
names(out_tab) <- c("Sample","Covariates","N (VPT)","N (FT)",
                    "Beta (Group)","t","p","Partial eta2","Cohen's d","R2")
writeData(wb, 1, out_tab, headerStyle = hs)
addStyle(wb, 1, cs, rows = 2:(nrow(out_tab)+1), cols = 1:ncol(out_tab), gridExpand = TRUE)
setColWidths(wb, 1, cols = 1:ncol(out_tab), widths = "auto")

# Sheet 2: Descriptives
addWorksheet(wb, "Descriptives")
desc <- do.call(rbind, lapply(
  list(c("Density","density"), c("Age at MRI","age5y"),
       c("Relative motion","FD"), c("eTIV","eTIV")),
  function(x) {
    do.call(rbind, lapply(c("VPT","FT"), function(g) {
      v <- df[[x[2]]][df$Group == g]
      data.frame(Variable = paste0(x[1], " (", g, ")"),
                 N = length(v), Mean = round(mean(v, na.rm = TRUE), 4),
                 SD = round(sd(v, na.rm = TRUE), 4),
                 Min = round(min(v, na.rm = TRUE), 4),
                 Max = round(max(v, na.rm = TRUE), 4))
    }))
  }))
writeData(wb, 2, desc, headerStyle = hs)
setColWidths(wb, 2, cols = 1:ncol(desc), widths = "auto")

# Sheet 3: Outliers
addWorksheet(wb, "Outliers")
out_info <- df[df$outlier, c("ID","Group","density","FD","eTIV","age5y")]
names(out_info) <- c("ID","Group","Density","Relative motion","eTIV","Age at MRI")
writeData(wb, 3, out_info, headerStyle = hs)
setColWidths(wb, 3, cols = 1:ncol(out_info), widths = "auto")

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat("Wrote Excel: ", OUT_XLSX, "\n", sep = "")

# ============================================================
# Figure 1: paneled boxplot + forest plot
# ============================================================

# Panel A: Boxplot + jitter (outliers in red)
df$pt_color <- ifelse(df$outlier, "Outlier (>2 SD)", "Included")

pA <- ggplot(df, aes(x = Group, y = density)) +
  geom_boxplot(aes(fill = Group), alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = pt_color, size = pt_color), width = 0.15, alpha = 0.6) +
  scale_fill_manual(values = c("FT" = "#2980b9", "VPT" = "#c0392b"), guide = "none") +
  scale_color_manual(values = c("Included" = "grey40", "Outlier (>2 SD)" = "#e31a1c"),
                     name = NULL) +
  scale_size_manual(values = c("Included" = 1.5, "Outlier (>2 SD)" = 3),
                    name = NULL) +
  annotate("segment", x = 1, xend = 2,
           y = max(df$density) + 0.008, yend = max(df$density) + 0.008,
           linewidth = 0.5) +
  annotate("text", x = 1.5, y = max(df$density) + 0.014,
           label = "p < 0.001", size = 3.5) +
  coord_cartesian(ylim = c(min(df$density) - 0.01, max(df$density) + 0.02)) +
  labs(x = NULL, y = "Connectome Density (unthresholded)",
       title = "A.  Group Comparison") +
  theme_classic(base_size = 12) +
  theme(legend.position.inside = c(0.15, 0.15),
        legend.position   = "inside",
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
        legend.text        = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 11))

# Panel B: Forest plot -- unadjusted and fully adjusted only
specs_forest <- list(
  "Unadjusted"     = density ~ Group,
  "Fully adjusted" = density ~ Group + age5y + sex + FD + eTIV
)

forest_df <- do.call(rbind, lapply(
  list(c("Full sample", "df"), c("Excl. outliers", "df_clean")),
  function(info) {
    ds <- get(info[2])
    do.call(rbind, lapply(names(specs_forest), function(nm) {
      m   <- lm(specs_forest[[nm]], data = ds)
      idx <- grep("Group", names(coef(m)))
      ci  <- confint(m)[idx, ]
      data.frame(Sample = info[1], Model = nm,
                 beta = coef(m)[idx], ci_lo = ci[1], ci_hi = ci[2],
                 p = summary(m)$coefficients[idx, "Pr(>|t|)"])
    }))
  }))

forest_df$Model  <- factor(forest_df$Model, levels = rev(names(specs_forest)))
forest_df$Sample <- factor(forest_df$Sample, levels = c("Full sample", "Excl. outliers"))
forest_df$p_lab  <- ifelse(forest_df$p < 0.001, "p < .001",
                           sprintf("p = %.3f", forest_df$p))

forest_df$y_num   <- as.numeric(forest_df$Model)
forest_df$y_dodge <- ifelse(forest_df$Sample == "Full sample",
                            forest_df$y_num - 0.15,
                            forest_df$y_num + 0.15)

pB <- ggplot(forest_df, aes(y = y_dodge, color = Sample)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = y_dodge), linewidth = 1) +
  geom_point(aes(x = beta), size = 3) +
  geom_text(aes(x = 0.002, label = p_lab),
            size = 3.2, hjust = 0, show.legend = FALSE, color = "grey25") +
  scale_y_continuous(breaks = 1:length(levels(forest_df$Model)),
                     labels = levels(forest_df$Model),
                     expand = expansion(mult = 0.3)) +
  scale_color_manual(values = c("Full sample"    = "#c0392b",
                                "Excl. outliers" = "#e67e22")) +
  coord_cartesian(xlim = c(min(forest_df$ci_lo) - 0.002,
                           max(forest_df$ci_hi) + 0.016)) +
  labs(x = expression(beta ~ "(VPT – FT)"), y = NULL, color = NULL,
       title = "B.  Sensitivity: Unadjusted vs. Fully Adjusted") +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        plot.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10))

# Combine panels
fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(1.2, 1), align = "v")
ggsave(OUT_FIG, fig, width = 8, height = 8, dpi = 300)

cat("Wrote figure: ", OUT_FIG, "\n", sep = "")
cat("\nDone.\n")