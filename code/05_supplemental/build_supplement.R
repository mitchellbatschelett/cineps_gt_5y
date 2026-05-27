################################################################################
#
#   CINEPS GT MANUSCRIPT --- SUPPLEMENT BUILDER
#
#   Generates Supplementary Figures and Tables from current RData files and
#   analytic data, using the same path, helper, and naming conventions as
#   `build_main_results.R`.
#
#.
#
################################################################################


# ============================================================================
# SECTION 0  CONFIGURATION
# ============================================================================

# --- Resolve repo root ----------------------------------------------------
# Tries three strategies, in order:
#   1. If this script was source()'d, derive its location from sys.frame
#      and step up to the repo root (two levels above code/05_supplemental).
#   2. If getwd() is already inside the repo, walk upward until we find one
#      with `data/analysis_ready/` and `code/`.
#   3. Fall back to hardcoded default; user can override before sourcing.
#
# Override by setting `repo_root` in the global environment BEFORE source():
#   repo_root <- "/full/path/to/repository"
#   source(".../build_supplement.R")

.is_repo_root <- function(p) {
  dir.exists(file.path(p, "data", "analysis_ready")) &&
    dir.exists(file.path(p, "code"))
}

.find_repo_root_from_path <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  for (i in 1:6) {
    if (.is_repo_root(p)) return(p)
    parent <- dirname(p)
    if (parent == p) break
    p <- parent
  }
  NULL
}

# Strategy 0: honor a pre-set repo_root from the global env
if (exists("repo_root", inherits = TRUE) &&
    is.character(repo_root) && length(repo_root) == 1 &&
    .is_repo_root(normalizePath(repo_root, winslash = "/", mustWork = FALSE))) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
} else {
  repo_root <- NULL
  
  # Strategy 1: derive from source()'d script location
  this_script <- tryCatch({
    sf <- sys.frames()
    paths <- character(0)
    for (fr in sf) {
      ofile <- tryCatch(get("ofile", envir = fr, inherits = FALSE),
                        error = function(e) NULL)
      if (!is.null(ofile) && is.character(ofile) && nzchar(ofile)) {
        paths <- c(paths, ofile)
      }
    }
    if (length(paths) > 0) paths[length(paths)] else NULL
  }, error = function(e) NULL)
  
  if (!is.null(this_script)) {
    repo_root <- .find_repo_root_from_path(dirname(this_script))
  }
  
  # Strategy 2: walk up from cwd
  if (is.null(repo_root)) {
    repo_root <- .find_repo_root_from_path(getwd())
  }
  
  if (is.null(repo_root)) {
    stop(
      "Could not locate repo root. Set it manually before source()-ing:\n",
      "  repo_root <- \"/full/path/to/repository\"\n",
      "  source(\".../build_supplement.R\")"
    )
  }
}
cat(sprintf("Repo root: %s\n", repo_root))

# --- Input paths -----------------------------------------------------------
data_xlsx  <- file.path(repo_root, "data", "analysis_ready",
                        "cohort_171VPT_45FT_postVQC.xlsx")
pilot_xlsx <- file.path(repo_root, "data", "intermediate",
                        "pilot_1000null_merged.xlsx")
frag_csv   <- file.path(repo_root, "data", "intermediate",
                        "fragmentation_by_density.csv")

# --- Output paths ----------------------------------------------------------
tables_dir  <- file.path(repo_root, "tables",  "supplement")
figures_dir <- file.path(repo_root, "figures", "supplement")
for (d in c(tables_dir, figures_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# --- Analytic density range (matches manuscript) --------------------------
density_min <- 11
density_max <- 100

# --- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(irr)
  library(patchwork)
})

# --- Plot theme  ------------------
theme_supp <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title         = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle      = element_text(size = base_size - 1, color = "grey30"),
      axis.title         = element_text(size = base_size),
      axis.text          = element_text(size = base_size - 2, color = "black"),
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 2),
      legend.position    = "right",
      legend.key.width   = unit(1.0, "cm"),
      strip.text         = element_text(size = base_size, face = "bold"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin        = margin(8, 12, 8, 8)
    )
}

# --- Color palette for metric pairs (consistent across SuppFigs 1 & 2) ----
pair_colors <- c(
  "ACC vs GE"        = "#1F77B4",  # blue
  "ACC vs Strength"  = "#D62728",  # red
  "GE vs Strength"   = "#2CA02C"   # green
)


# ============================================================================
# SECTION 1  HELPERS
# ============================================================================

# Save a figure as both PNG (300 dpi) and PDF (cairo) into figures/supplement/
save_fig <- function(plot_obj, name, width_in, height_in) {
  png_path <- file.path(figures_dir, paste0(name, ".png"))
  pdf_path <- file.path(figures_dir, paste0(name, ".pdf"))
  ggsave(png_path, plot = plot_obj, width = width_in, height = height_in,
         dpi = 300, units = "in", bg = "white")
  ggsave(pdf_path, plot = plot_obj, width = width_in, height = height_in,
         units = "in", device = cairo_pdf)
  cat(sprintf("  Saved: %s.png + %s.pdf\n", name, name))
  invisible(plot_obj)
}

# Identify density columns matching a metric prefix within [d_min, d_max].
get_density_cols <- function(df, metric_prefix, d_min = density_min, d_max = density_max) {
  pat <- paste0("^", metric_prefix, "_[0-9]+(\\.[0-9]+)?$")
  candidates <- grep(pat, names(df), value = TRUE)
  d_vals <- as.numeric(sub(paste0("^", metric_prefix, "_"), "", candidates))
  keep <- !is.na(d_vals) & d_vals >= d_min & d_vals <= d_max
  candidates <- candidates[keep]
  d_vals     <- d_vals[keep]
  ord <- order(d_vals)
  list(cols = candidates[ord], densities = d_vals[ord])
}

# Pearson r at every density for two metric prefixes. Returns a tibble:
# density, r, n  (n = pairwise complete observations).
pairwise_r_by_density <- function(df, prefix_x, prefix_y,
                                  d_min = density_min, d_max = density_max) {
  ax <- get_density_cols(df, prefix_x, d_min, d_max)
  ay <- get_density_cols(df, prefix_y, d_min, d_max)
  
  shared <- intersect(ax$densities, ay$densities)
  if (length(shared) == 0) {
    stop(sprintf("No overlapping densities between %s and %s in range [%g, %g].",
                 prefix_x, prefix_y, d_min, d_max))
  }
  
  rs <- numeric(length(shared))
  ns <- integer(length(shared))
  for (i in seq_along(shared)) {
    d  <- shared[i]
    cx <- ax$cols[match(d, ax$densities)]
    cy <- ay$cols[match(d, ay$densities)]
    x  <- as.numeric(df[[cx]])
    y  <- as.numeric(df[[cy]])
    ok <- complete.cases(x, y)
    ns[i] <- sum(ok)
    rs[i] <- if (ns[i] >= 3) suppressWarnings(cor(x[ok], y[ok], method = "pearson")) else NA_real_
  }
  
  tibble(density = shared, r = rs, n = ns)
}


# ============================================================================
# SECTION 2  DATA LOADING
# ============================================================================

cat("\n========================================\n")
cat("Loading analytic xlsx\n")
cat("========================================\n\n")

if (!file.exists(data_xlsx)) {
  stop(sprintf("Analytic xlsx not found at:\n  %s", data_xlsx))
}

df <- as.data.frame(read_excel(data_xlsx))
cat(sprintf("Loaded: %d rows x %d cols\n", nrow(df), ncol(df)))

if ("Group" %in% names(df)) {
  df$Group <- factor(df$Group)
  grp_tab  <- table(df$Group)
  cat(sprintf("Group counts: %s (Group=1 should be VPT, Group=0 should be FT)\n",
              paste(sprintf("%s=%d", names(grp_tab), as.integer(grp_tab)),
                    collapse = ", ")))
}


# ============================================================================
# SECTION N1  NOTE 1 --- METRIC NORMALIZATION DEMONSTRATION
# ----------------------------------------------------------------------------
# Generates:
#   Supplementary Figure 1  --- Pre-normalization r(density), 3 metric pairs
#   Supplementary Figure 2  --- Post-normalization r(density), same 3 pairs
#   Supplementary Table 1   --- Mean / min / max r per pair, per stage
#
# Method:
#   At each density in [11, 100]% (steps of 0.25, so 357 grid points),
#   compute pairwise Pearson r between metric columns across the n=216
#   participants. Three pairs are reported:
#       (1) ACC vs GE
#       (2) ACC vs Strength
#       (3) GE  vs Strength
#   Pre-normalization uses the raw weighted metrics: ACC, GE, str.
#   Post-normalization uses the random-network-normalized counterparts:
#   rand_norm_wei_ACC, rand_norm_wei_GE, str (strength is unchanged by
#   null-model normalization, so it appears in both stages).
# ============================================================================

cat("\n========================================\n")
cat("Note 1: Normalization Demonstration\n")
cat("========================================\n\n")

# --- Compute pairwise r curves --------------------------------------------
cat("Computing pre-normalization r(density)...\n")
pre_acc_ge <- pairwise_r_by_density(df, "ACC", "GE")    %>% mutate(pair = "ACC vs GE")
pre_acc_st <- pairwise_r_by_density(df, "ACC", "str")   %>% mutate(pair = "ACC vs Strength")
pre_ge_st  <- pairwise_r_by_density(df, "GE",  "str")   %>% mutate(pair = "GE vs Strength")
pre_df <- bind_rows(pre_acc_ge, pre_acc_st, pre_ge_st) %>%
  mutate(pair = factor(pair, levels = names(pair_colors)))

cat("Computing post-normalization r(density)...\n")
post_acc_ge <- pairwise_r_by_density(df, "rand_norm_wei_ACC", "rand_norm_wei_GE")  %>% mutate(pair = "ACC vs GE")
post_acc_st <- pairwise_r_by_density(df, "rand_norm_wei_ACC", "str")               %>% mutate(pair = "ACC vs Strength")
post_ge_st  <- pairwise_r_by_density(df, "rand_norm_wei_GE",  "str")               %>% mutate(pair = "GE vs Strength")
post_df <- bind_rows(post_acc_ge, post_acc_st, post_ge_st) %>%
  mutate(pair = factor(pair, levels = names(pair_colors)))

cat(sprintf("  Pre-norm: %d density points per pair (n per density = %d)\n",
            nrow(pre_acc_ge), unique(pre_acc_ge$n)[1]))
cat(sprintf("  Post-norm: %d density points per pair\n\n",
            nrow(post_acc_ge)))


# --- Supplementary Table 1: pre/post r summary ----------------------------
note1_summary <- bind_rows(
  pre_df  %>% group_by(pair) %>%
    summarise(stage = "Pre-normalization",
              mean_r = mean(r, na.rm = TRUE),
              min_r  = min(r,  na.rm = TRUE),
              max_r  = max(r,  na.rm = TRUE),
              .groups = "drop"),
  post_df %>% group_by(pair) %>%
    summarise(stage = "Post-normalization",
              mean_r = mean(r, na.rm = TRUE),
              min_r  = min(r,  na.rm = TRUE),
              max_r  = max(r,  na.rm = TRUE),
              .groups = "drop")
) %>%
  mutate(
    stage = factor(stage, levels = c("Pre-normalization", "Post-normalization"))
  ) %>%
  arrange(stage, pair) %>%
  transmute(
    Stage      = as.character(stage),
    `Metric pair` = as.character(pair),
    `Mean r`   = sprintf("%.3f", mean_r),
    `Min r`    = sprintf("%.3f", min_r),
    `Max r`    = sprintf("%.3f", max_r)
  )

cat("Supplementary Table 1: pairwise correlation summary\n\n")
print(note1_summary)

write.csv(note1_summary,
          file.path(tables_dir, "SuppTbl1_normalization_r_summary.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path("tables", "supplement",
                      "SuppTbl1_normalization_r_summary.csv")))


# --- Shared y-axis range so Figs 1 and 2 are visually comparable ----------
y_min <- floor(  min(c(pre_df$r, post_df$r), na.rm = TRUE) * 10) / 10
y_max <- ceiling(max(c(pre_df$r, post_df$r), na.rm = TRUE) * 10) / 10
y_breaks <- seq(y_min, y_max, by = 0.25)


# --- Plot builder (Figs 1 and 2 differ only in input data) ----------------
make_norm_fig <- function(plot_df, title) {
  ggplot(plot_df, aes(x = density, y = r, color = pair)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pair_colors, name = "Metric pair") +
    scale_x_continuous(breaks = seq(20, 100, by = 20),
                       limits = c(density_min, density_max),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = y_breaks,
                       limits = c(y_min, y_max),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Network density (%)",
         y = "Pearson r",
         title = title) +
    theme_supp()
}

p_fig1 <- make_norm_fig(
  pre_df,
  "Supplementary Figure 1. Pre-normalization metric correlations across density"
)
p_fig2 <- make_norm_fig(
  post_df,
  "Supplementary Figure 2. Post-normalization metric correlations across density"
)

cat("\nWriting figures...\n")
save_fig(p_fig1, "SuppFig1_pre_normalization_r_by_density",  width_in = 7.0, height_in = 4.0)
save_fig(p_fig2, "SuppFig2_post_normalization_r_by_density", width_in = 7.0, height_in = 4.0)


cat("\n========================================\n")
cat("Note 1: DONE\n")
cat("========================================\n")


# ============================================================================
# SECTION N2  NOTE 2 --- NULL MODEL VALIDATION
# ----------------------------------------------------------------------------
# Generates:
#   Supplementary Figure 3  --- 3-panel Pearson r curves comparing 100-null
#                                vs 1000-null normalization across density,
#                                one panel per metric (GE, ACC, SW)
#   Supplementary Table 2   --- Per-metric mean / min / max of Pearson r AND
#                                ICC(2,1) across the analytic density range
#
# Method:
#   The pilot xlsx contains 1000-null estimates for 30 randomly selected
#   subjects (15 VPT + 15 FT, seed=42), held out by ID from the main 100-null
#   analytic xlsx. At every density in [11, 100]% we compute Pearson r and
#   ICC(2,1) between the 100-null and 1000-null estimates across the 30
#   subjects, then summarize across density. Strength is unchanged by
#   null-model normalization and is therefore not reported here. L is
#   dropped per manuscript scope.
# ============================================================================

cat("\n========================================\n")
cat("Note 2: Null Model Validation\n")
cat("========================================\n\n")

if (!file.exists(pilot_xlsx)) {
  stop(sprintf("Pilot 1000-null xlsx not found at:\n  %s", pilot_xlsx))
}

cat("Loading pilot 1000-null xlsx...\n")
df_pilot <- as.data.frame(read_excel(pilot_xlsx))

# Strip "sub-" prefix if present (the current merge_pilot_validation.m strips
# it already; this is a belt-and-suspenders pass in case an older pilot xlsx
# is used).
df_pilot$ID <- gsub("^sub-", "", as.character(df_pilot$ID))

pilot_n   <- nrow(df_pilot)
pilot_vpt <- sum(df_pilot$Group == 1, na.rm = TRUE)
pilot_ft  <- sum(df_pilot$Group == 0, na.rm = TRUE)
cat(sprintf("  Pilot N: %d (VPT=%d, FT=%d)\n", pilot_n, pilot_vpt, pilot_ft))

# --- Align main xlsx subjects to pilot ID order ---------------------------
df$ID <- as.character(df$ID)
pilot_ids <- df_pilot$ID
df_main_sub <- df[match(pilot_ids, df$ID), , drop = FALSE]

n_matched <- sum(!is.na(df_main_sub$ID))
if (n_matched < pilot_n) {
  warning(sprintf("Note 2: %d pilot IDs missing from main xlsx", pilot_n - n_matched))
}
cat(sprintf("  Matched in main xlsx: %d / %d\n\n", n_matched, pilot_n))


# --- Metrics to validate (matches manuscript: GE, ACC, SW) -----------------
n2_metrics <- c("rand_norm_wei_GE", "rand_norm_wei_ACC", "rand_norm_wei_SW")
n2_metric_labels <- c(
  "rand_norm_wei_GE"  = "Global Efficiency (GE)",
  "rand_norm_wei_ACC" = "Average Clustering (ACC)",
  "rand_norm_wei_SW"  = "Small-Worldness (SW)"
)

# --- Per-density r AND ICC(2,1) for each metric ---------------------------
n2_density_list <- list()

for (m in n2_metrics) {
  ml <- n2_metric_labels[m]
  cat(sprintf("--- %s ---\n", ml))
  
  ax_main  <- get_density_cols(df_main_sub, m, density_min, density_max)
  ax_pilot <- get_density_cols(df_pilot,    m, density_min, density_max)
  
  shared_d <- intersect(ax_main$densities, ax_pilot$densities)
  if (length(shared_d) == 0) {
    warning(sprintf("Note 2: no shared densities for %s", m))
    next
  }
  
  # Get aligned column lists in shared-density order
  main_cols  <- ax_main$cols [match(shared_d, ax_main$densities)]
  pilot_cols <- ax_pilot$cols[match(shared_d, ax_pilot$densities)]
  
  main_mat  <- as.matrix(df_main_sub[, main_cols,  drop = FALSE])
  pilot_mat <- as.matrix(df_pilot[,    pilot_cols, drop = FALSE])
  mode(main_mat)  <- "numeric"
  mode(pilot_mat) <- "numeric"
  
  rvals   <- numeric(length(shared_d))
  iccvals <- numeric(length(shared_d))
  
  for (j in seq_along(shared_d)) {
    o_d <- main_mat[,  j]
    p_d <- pilot_mat[, j]
    v <- is.finite(o_d) & is.finite(p_d)
    if (sum(v) >= 3 && sd(o_d[v]) > 0 && sd(p_d[v]) > 0) {
      rvals[j] <- cor(o_d[v], p_d[v])
      icc_obj <- suppressWarnings(
        icc(data.frame(m1 = o_d[v], m2 = p_d[v]),
            model = "twoway", type = "agreement", unit = "single")
      )
      iccvals[j] <- icc_obj$value
    } else {
      rvals[j]   <- NA_real_
      iccvals[j] <- NA_real_
    }
  }
  
  n2_density_list[[m]] <- tibble(
    metric  = ml,
    density = shared_d,
    r       = rvals,
    ICC     = iccvals
  )
  
  cat(sprintf("  r: mean=%.4f, min=%.4f, max=%.4f\n",
              mean(rvals, na.rm = TRUE), min(rvals, na.rm = TRUE), max(rvals, na.rm = TRUE)))
  cat(sprintf("  ICC: mean=%.4f, min=%.4f, max=%.4f\n",
              mean(iccvals, na.rm = TRUE), min(iccvals, na.rm = TRUE), max(iccvals, na.rm = TRUE)))
}

n2_density_df <- bind_rows(n2_density_list) %>%
  mutate(metric = factor(metric, levels = unname(n2_metric_labels)))


# --- Supplementary Table 2: per-metric summary across density --------------
note2_summary <- n2_density_df %>%
  group_by(metric) %>%
  summarise(
    mean_r   = mean(r,   na.rm = TRUE),
    min_r    = min(r,    na.rm = TRUE),
    max_r    = max(r,    na.rm = TRUE),
    mean_ICC = mean(ICC, na.rm = TRUE),
    min_ICC  = min(ICC,  na.rm = TRUE),
    max_ICC  = max(ICC,  na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  transmute(
    Metric     = as.character(metric),
    `Mean r`   = sprintf("%.3f", mean_r),
    `Min r`    = sprintf("%.3f", min_r),
    `Max r`    = sprintf("%.3f", max_r),
    `Mean ICC` = sprintf("%.3f", mean_ICC),
    `Min ICC`  = sprintf("%.3f", min_ICC),
    `Max ICC`  = sprintf("%.3f", max_ICC)
  )

cat("\nSupplementary Table 2: null model convergence summary\n\n")
print(note2_summary)

write.csv(note2_summary,
          file.path(tables_dir, "SuppTbl2_null_validation_summary.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path("tables", "supplement",
                      "SuppTbl2_null_validation_summary.csv")))


# --- Supplementary Figure 3: 3-panel r-by-density curves ------------------
# Y-axis: low end determined by data, capped at most 0.95 to keep visual context
y_lo_n2 <- min(0.95, floor(min(n2_density_df$r, na.rm = TRUE) * 100) / 100)

p_fig3 <- ggplot(n2_density_df, aes(x = density, y = r)) +
  geom_hline(yintercept = 1.00, color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0.99, color = "#D62728", linetype = "dashed", linewidth = 0.4) +
  geom_line(color = "#1F77B4", linewidth = 0.9) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20),
                     limits = c(density_min, density_max),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(y_lo_n2, 1.005),
                     breaks = seq(round(y_lo_n2, 2), 1.00, by = 0.01),
                     expand = expansion(mult = c(0.01, 0.005))) +
  labs(x = "Network density (%)",
       y = "Pearson r (across 30 subjects)",
       title = "Supplementary Figure 3. 100-null vs 1000-null convergence by density") +
  theme_supp()

cat("\nWriting figures...\n")
save_fig(p_fig3, "SuppFig3_null_validation_r_by_density",
         width_in = 9.0, height_in = 3.5)


cat("\n========================================\n")
cat("Note 2: DONE\n")
cat("========================================\n")


# ============================================================================
# SECTION N3  NOTE 3 --- DENSITY LOWER BOUND DETERMINATION
# ----------------------------------------------------------------------------
# Underlying fragmentation analysis was performed in MATLAB
# (network_fragmentation.m, in code/05_supplemental/). This section reads the
# MATLAB-generated CSV and renders Supplementary Figure 4.
#
# Generates:
#   Supplementary Figure 4 --- Percentage of subjects with single connected
#                              component as a function of network density.
# ============================================================================

cat("\n========================================\n")
cat("Note 3: Density Lower Bound Determination\n")
cat("========================================\n\n")

if (!file.exists(frag_csv)) {
  stop(sprintf("Fragmentation CSV not found at:\n  %s", frag_csv))
}

# --- Load MATLAB-generated fragmentation data -----------------------------
cat("Loading fragmentation_by_density.csv...\n")
frag_df <- read.csv(frag_csv, stringsAsFactors = FALSE)
cat(sprintf("  Loaded: %d densities, %d columns\n", nrow(frag_df), ncol(frag_df)))

# Identify the density at which 100% of subjects have a single connected component
threshold_idx <- which(frag_df$pct_single_comp_all >= 100)[1]
threshold_density <- if (!is.na(threshold_idx)) frag_df$density_pct[threshold_idx] else NA
cat(sprintf("  All-connected threshold (all subjects): %.2f%%\n\n", threshold_density))

# --- Reshape for ggplot: one line per group (All / VPT / FT) --------------
n3_df <- data.frame(
  density = rep(frag_df$density_pct, 3),
  pct     = c(frag_df$pct_single_comp_all,
              frag_df$pct_single_comp_vpt,
              frag_df$pct_single_comp_ft),
  group   = factor(rep(c("All", "VPT", "FT"), each = nrow(frag_df)),
                   levels = c("All", "VPT", "FT"))
)

# Subset x-range for visualization (curve plateaus well before 25%; show 0-25%)
xmax_n3 <- 25

# Color palette: keep neutral for "All", same VPT/FT colors used elsewhere
n3_colors <- c(
  "All" = "black",
  "VPT" = "#1F77B4",
  "FT"  = "#D62728"
)

# --- Build figure ---------------------------------------------------------
p_fig4 <- ggplot(subset(n3_df, density <= xmax_n3),
                 aes(x = density, y = pct, color = group, linetype = group)) +
  geom_vline(xintercept = threshold_density, color = "grey50",
             linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 100, color = "grey60",
             linetype = "dotted", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  annotate("text",
           x = threshold_density + 0.4, y = 8,
           label = sprintf("All connected: %.1f%%", threshold_density),
           color = "grey25", size = 3.2, hjust = 0, fontface = "italic") +
  scale_color_manual(values = n3_colors, name = "Group") +
  scale_linetype_manual(values = c("All" = "solid", "VPT" = "solid", "FT" = "solid"),
                        guide = "none") +
  scale_x_continuous(breaks = seq(0, xmax_n3, by = 5),
                     limits = c(0, xmax_n3),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     limits = c(0, 102),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Network density (%)",
       y = "% of subjects with single connected component",
       title = "Supplementary Figure 4. Network fragmentation across density thresholds") +
  theme_supp()

cat("Writing figures...\n")
save_fig(p_fig4, "SuppFig4_fragmentation_pct_single_component",
         width_in = 7.0, height_in = 4.0)


cat("\n========================================\n")
cat("Note 3: DONE\n")
cat("========================================\n")


# ============================================================================
# SECTION N4  NOTE 4 --- GROUP DIFFERENCE MODEL DETAILS
# ----------------------------------------------------------------------------
# Generates (across sub-sections):
#   4a  Supplementary Figure 5  --- VPT vs FT raw mean +/- SD curves (6 panels)
#   4b  Supplementary Table 3   --- Density single-cov ANCOVA (8 rows: full +
#                                    outlier-excluded samples)
#   4c  Supplementary Table 4   --- Per-IV PFFR output (84 rows)  [pending]
#   4d  Supplementary Table 5   --- K-basis diagnostics (36 rows) [pending]
#   4e  Supplementary Figs 6-10 --- Covariate beta(d) curves      [pending]
#
# 4a and 4b read from the analytic xlsx already loaded in Section 2.
# 4c reads predictor_tests.csv + summary.csv from results/fda_group_differences/.
# 4d-4e load fully-adjusted PFFR RData from the same directories.
# ============================================================================

cat("\n========================================\n")
cat("Note 4: Group Difference Model Details\n")
cat("========================================\n\n")

# Per-metric subject exclusions
n4_metric_excl <- list(
  "str"               = c("128", "713"),
  "rand_norm_wei_GE"  = c("659"),
  "rand_norm_wei_ACC" = character(0),
  "rand_norm_wei_SW"  = c("309", "321", "8155"),
  "GE"                = c("128", "713"),
  "ACC"               = c("128", "713")
)

# --- Metric specification: column prefix in xlsx + display label + panel order ---
n4_metrics <- list(
  list(prefix = "str",                label = "Strength"),
  list(prefix = "rand_norm_wei_GE",   label = "Normalized GE"),
  list(prefix = "rand_norm_wei_ACC",  label = "Normalized ACC"),
  list(prefix = "rand_norm_wei_SW",   label = "Small-Worldness"),
  list(prefix = "GE",                 label = "Raw GE"),
  list(prefix = "ACC",                label = "Raw ACC")
)
n4_panel_levels <- vapply(n4_metrics, function(x) x$label, character(1))


# ----------------------------------------------------------------------------
# 4a  Supplementary Figure 5 --- VPT vs FT mean +/- SD curves
# ----------------------------------------------------------------------------

cat("--- 4a: Building Supplementary Figure 5 ---\n\n")

n4_curves <- list()

for (mspec in n4_metrics) {
  ax <- get_density_cols(df, mspec$prefix, density_min, density_max)
  if (length(ax$cols) == 0) {
    cat(sprintf("  [skip] %s: no density columns in [%g, %g]%%\n",
                mspec$label, density_min, density_max))
    next
  }
  
  # Apply metric-specific exclusions
  excl_ids <- n4_metric_excl[[mspec$prefix]]
  if (length(excl_ids) > 0) {
    keep_rows <- !(as.character(df$ID) %in% as.character(excl_ids))
  } else {
    keep_rows <- rep(TRUE, nrow(df))
  }
  
  df_m  <- df[keep_rows, , drop = FALSE]
  mat   <- as.matrix(df_m[, ax$cols, drop = FALSE])
  group <- df_m$Group
  
  # Per-metric complete-cases: keep subjects with non-NA at every density column
  ok <- complete.cases(mat)
  mat_ok   <- mat[ok, , drop = FALSE]
  group_ok <- group[ok]
  
  # Group is coded 1 = VPT, 0 = FT
  is_vpt <- group_ok == 1
  is_ft  <- group_ok == 0
  
  vpt_mean <- colMeans(mat_ok[is_vpt, , drop = FALSE])
  vpt_sd   <- apply(mat_ok[is_vpt, , drop = FALSE], 2, sd)
  ft_mean  <- colMeans(mat_ok[is_ft,  , drop = FALSE])
  ft_sd    <- apply(mat_ok[is_ft,  , drop = FALSE], 2, sd)
  
  n4_curves[[mspec$label]] <- bind_rows(
    tibble(metric = mspec$label, density = ax$densities,
           mean = vpt_mean, sd = vpt_sd, group = "VPT", n = sum(is_vpt)),
    tibble(metric = mspec$label, density = ax$densities,
           mean = ft_mean,  sd = ft_sd,  group = "FT",  n = sum(is_ft))
  )
  
  cat(sprintf("  %-18s n_total=%d (VPT=%d, FT=%d)  [excluded %d IDs]\n",
              mspec$label, sum(ok), sum(is_vpt), sum(is_ft), length(excl_ids)))
}

n4_curves_df <- bind_rows(n4_curves) %>%
  mutate(
    metric = factor(metric, levels = n4_panel_levels),
    group  = factor(group,  levels = c("VPT", "FT"))
  )

# 6 panels (2 rows x 3 cols), VPT vs FT mean lines + 1-SD ribbons
group_colors <- c("VPT" = "#1F77B4", "FT" = "#D62728")

p_fig5 <- ggplot(n4_curves_df,
                 aes(x = density, y = mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ metric, ncol = 3, scales = "free_y") +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_fill_manual( values = group_colors, guide = "none") +
  scale_x_continuous(breaks = seq(20, 100, by = 20),
                     limits = c(density_min, density_max),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Network density (%)",
       y = "Metric value",
       title = "Supplementary Figure 5. VPT vs FT metric-by-density curves") +
  theme_supp() +
  theme(legend.position = "top")

cat("\nWriting figure...\n")
save_fig(p_fig5, "SuppFig5_raw_curves_VPT_vs_FT",
         width_in = 9.0, height_in = 6.0)


# ----------------------------------------------------------------------------
# 4b  Supplementary Table 3 --- Density single-covariate ANCOVA
# ----------------------------------------------------------------------------
# Full sample (n=171/45) plus outlier-excluded sample (n=165/45). Four
# single-covariate ANCOVA fits per sample: Group + {Age at MRI | Sex |
# Relative motion | TIV}. Models predict den_100.00 (the unthresholded
# matrix density).
#
# Reports per model: Beta(Group), t(Group), p(Group), partial eta-squared
# for Group, and overall R^2. The outlier-excluded sample drops VPT subjects
# whose density is more than 2 SD below the VPT mean.
# ----------------------------------------------------------------------------

cat("\n--- 4b: Building Supplementary Table 3 ---\n\n")

# Subset to complete cases on all variables used
n4_cov_vars <- c("Group", "den_100.00", "age_at_5y_mri", "sex", "Rel_Motion", "eTIV")
miss_any <- !apply(df[, n4_cov_vars], 1, function(r) all(!is.na(r)))
df_t3 <- df[!miss_any, , drop = FALSE]
df_t3$Group <- as.numeric(as.character(df_t3$Group))   # 1=VPT, 0=FT, numeric for OLS

cat(sprintf("  Sup Tbl 3 input n=%d (VPT=%d, FT=%d)\n",
            nrow(df_t3), sum(df_t3$Group == 1), sum(df_t3$Group == 0)))

# Outlier exclusion: VPT subjects with density > 2 SD below VPT mean
vpt_density <- df_t3$`den_100.00`[df_t3$Group == 1]
vpt_mean_d  <- mean(vpt_density)
vpt_sd_d    <- sd(vpt_density)
outlier_thr <- vpt_mean_d - 2 * vpt_sd_d
out_mask    <- df_t3$Group == 1 & df_t3$`den_100.00` < outlier_thr
df_t3_excl  <- df_t3[!out_mask, , drop = FALSE]

cat(sprintf("  VPT density: mean=%.4f, SD=%.4f, 2SD-below=%.4f, n outliers=%d\n",
            vpt_mean_d, vpt_sd_d, outlier_thr, sum(out_mask)))
cat(sprintf("  After exclusion: n=%d (VPT=%d, FT=%d)\n\n",
            nrow(df_t3_excl), sum(df_t3_excl$Group == 1), sum(df_t3_excl$Group == 0)))

# Helper: fit Group + one covariate ANCOVA on density, return formatted stats
fit_density_ancova <- function(data, covariate) {
  # Full model: y ~ Group + covariate
  fm <- as.formula(sprintf("`den_100.00` ~ Group + %s", covariate))
  full <- lm(fm, data = data)
  sm   <- summary(full)
  cf   <- sm$coefficients
  
  beta_g <- cf["Group", "Estimate"]
  t_g    <- cf["Group", "t value"]
  p_g    <- cf["Group", "Pr(>|t|)"]
  r2     <- sm$r.squared
  
  # Partial eta-squared for Group = SS_Group / (SS_Group + SS_resid_full)
  # SS_Group is the incremental SS from dropping Group.
  fm_red <- as.formula(sprintf("`den_100.00` ~ %s", covariate))
  red    <- lm(fm_red, data = data)
  ss_res_full <- sum(residuals(full)^2)
  ss_res_red  <- sum(residuals(red)^2)
  ss_group    <- ss_res_red - ss_res_full
  peta2       <- ss_group / (ss_group + ss_res_full)
  
  list(beta = beta_g, t = t_g, p = p_g, peta2 = peta2, r2 = r2)
}

# Covariate spec: (column name in xlsx, publication label)
n4_t3_covs <- list(
  list(var = "age_at_5y_mri", label = "Age at MRI"),
  list(var = "sex",           label = "Sex"),
  list(var = "Rel_Motion",    label = "Relative motion"),
  list(var = "eTIV",          label = "TIV")
)

fmt_p_t3 <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<.001")
  sprintf("%.3f", p)
}

build_t3_block <- function(data, sample_label, n_vpt, n_ft) {
  rows <- list()
  for (cv in n4_t3_covs) {
    r <- fit_density_ancova(data, cv$var)
    rows[[length(rows) + 1]] <- tibble(
      Sample                 = sample_label,
      `Covariate added`      = cv$label,
      `n (VPT/FT)`           = sprintf("%d/%d", n_vpt, n_ft),
      `Beta (VPT-FT)`        = sprintf("%.3f", r$beta),
      `t`                    = sprintf("%.2f", r$t),
      `p`                    = fmt_p_t3(r$p),
      `Partial eta-squared`  = sprintf("%.3f", r$peta2),
      `R-squared`            = sprintf("%.3f", r$r2)
    )
  }
  bind_rows(rows)
}

t3_full <- build_t3_block(df_t3,      "Full sample",
                          sum(df_t3$Group == 1), sum(df_t3$Group == 0))
t3_excl <- build_t3_block(df_t3_excl, "Excl. outliers (2 SD)",
                          sum(df_t3_excl$Group == 1), sum(df_t3_excl$Group == 0))

# Preserve covariate order from the supplement (Age, Sex, Motion, TIV)
n4_t3_cov_order <- vapply(n4_t3_covs, function(x) x$label, character(1))

note4_t3 <- bind_rows(t3_full, t3_excl) %>%
  mutate(
    Sample            = factor(Sample,
                               levels = c("Full sample", "Excl. outliers (2 SD)")),
    `Covariate added` = factor(`Covariate added`, levels = n4_t3_cov_order)
  ) %>%
  arrange(Sample, `Covariate added`)

cat("Supplementary Table 3: density single-covariate ANCOVAs\n\n")
print(note4_t3, n = Inf)

write.csv(note4_t3,
          file.path(tables_dir, "SuppTbl3_density_single_cov_ANCOVA.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path("tables", "supplement",
                      "SuppTbl3_density_single_cov_ANCOVA.csv")))


# ----------------------------------------------------------------------------
# 4c  Supplementary Table 4 --- Per-IV PFFR output across 6 covariate configs
# ----------------------------------------------------------------------------
# For each of 6 metrics x 6 covariate configs, reads the precomputed
# predictor_tests.csv (per-term test statistic, p, functional VE/dev-explained,
# semi-partial r) and summary.csv (family, N) from
# results/fda_group_differences/<config_dir>/.
#
# SW models use the scat (scaled t) family, so their Test_Statistic is chi^2
# (not F) and their "functional_ve" is deviance explained (not variance
# explained). The table footnote documents this; the values are reported on
# the same numerical scale (%, x100 from the [0,1] CSV value).
# ----------------------------------------------------------------------------

cat("\n--- 4c: Building Supplementary Table 4 ---\n\n")

# --- Config spec: directory suffix + display label, in supplement row order ---
n4_t4_configs <- list(
  list(suffix = "Group_11-100_fullsample",
       label  = "Unadjusted"),
  list(suffix = "Group_eTIV_11-100_fullsample",
       label  = "TIV"),
  list(suffix = "Group_age_at_5y_mri_11-100_fullsample",
       label  = "Age at MRI"),
  list(suffix = "Group_sex_11-100_fullsample",
       label  = "Sex"),
  list(suffix = "Group_Rel_Motion_11-100_fullsample",
       label  = "Relative motion"),
  list(suffix = "Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       label  = "Fully adjusted")
)
n4_t4_config_levels <- vapply(n4_t4_configs, function(x) x$label, character(1))

# --- Metric spec: directory prefix + display label, supplement row order ---
n4_t4_metrics <- list(
  list(prefix = "str",                label = "Strength"),
  list(prefix = "rand_norm_wei_GE",   label = "Normalized GE"),
  list(prefix = "rand_norm_wei_ACC",  label = "Normalized ACC"),
  list(prefix = "rand_norm_wei_SW",   label = "Small-Worldness"),
  list(prefix = "GE",                 label = "Raw GE"),
  list(prefix = "ACC",                label = "Raw ACC")
)
n4_t4_metric_levels <- vapply(n4_t4_metrics, function(x) x$label, character(1))

# --- Predictor display labels (matches supplement Table 4 conventions) ---
n4_t4_predictor_labels <- c(
  "Group"          = "Group (VPT vs FT)",
  "eTIV"           = "TIV",
  "age_at_5y_mri"  = "Age at MRI",
  "sex"            = "Sex",
  "Rel_Motion"     = "Relative motion"
)

# --- p-value formatter ---
fmt_p_t4 <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# --- Read one (metric, config) directory's predictor_tests + summary ---
read_t4_row_block <- function(metric_prefix, config_suffix, metric_label, config_label) {
  dir_name <- sprintf("%s_FDA_%s", metric_prefix, config_suffix)
  dir_path <- file.path(repo_root, "results", "fda_group_differences", dir_name)
  
  if (!dir.exists(dir_path)) {
    warning(sprintf("Missing directory: %s", dir_path))
    return(NULL)
  }
  
  pt_path  <- file.path(dir_path, sprintf("%s_predictor_tests.csv", metric_prefix))
  sum_path <- file.path(dir_path, sprintf("%s_summary.csv",         metric_prefix))
  
  if (!file.exists(pt_path)) {
    warning(sprintf("Missing predictor_tests.csv: %s", pt_path))
    return(NULL)
  }
  if (!file.exists(sum_path)) {
    warning(sprintf("Missing summary.csv: %s", sum_path))
    return(NULL)
  }
  
  pt   <- read.csv(pt_path,  stringsAsFactors = FALSE)
  smry <- read.csv(sum_path, stringsAsFactors = FALSE)
  
  # Family: "scat" or "scaled t" -> SW-style (chi^2 + deviance explained); else F
  fam <- tolower(as.character(smry$Family[1]))
  is_scat <- fam %in% c("scat", "scaled t")
  
  # Build display rows
  tibble(
    Metric                 = metric_label,
    `Covariate model`      = config_label,
    Predictor              = ifelse(pt$Predictor %in% names(n4_t4_predictor_labels),
                                    n4_t4_predictor_labels[pt$Predictor],
                                    pt$Predictor),
    is_scat                = is_scat,
    Test_Statistic_num     = pt$Test_Statistic,
    p_value                = pt$p_value,
    functional_ve_pct      = pt$functional_ve * 100,
    semi_partial_r         = pt$semi_partial_r,
    .rows = nrow(pt)
  )
}

# --- Loop over all metric x config combinations ---
t4_blocks <- list()
n_dirs_total   <- length(n4_t4_metrics) * length(n4_t4_configs)
n_dirs_found   <- 0

for (mspec in n4_t4_metrics) {
  for (cspec in n4_t4_configs) {
    rb <- read_t4_row_block(mspec$prefix, cspec$suffix, mspec$label, cspec$label)
    if (is.null(rb)) next
    t4_blocks[[length(t4_blocks) + 1]] <- rb
    n_dirs_found <- n_dirs_found + 1
  }
}

cat(sprintf("  Loaded %d / %d (metric x config) directories\n",
            n_dirs_found, n_dirs_total))

t4_raw <- bind_rows(t4_blocks)

# --- Format display ---
# Test statistic: SW gets superscript-b marker (chi^2 footnote); others as F.
note4_t4 <- t4_raw %>%
  mutate(
    `Test statistic` = ifelse(
      is_scat,
      paste0(sprintf("%.2f", Test_Statistic_num), "\u1d47"),   # superscript b
      sprintf("%.2f", Test_Statistic_num)
    ),
    `p`              = vapply(p_value, fmt_p_t4, character(1)),
    `Functional Variance Explained (%)` = sprintf("%.2f", functional_ve_pct),
    `Semi-partial r` = sprintf("%.3f", semi_partial_r),
    Metric            = factor(Metric,            levels = n4_t4_metric_levels),
    `Covariate model` = factor(`Covariate model`, levels = n4_t4_config_levels)
  ) %>%
  arrange(Metric, `Covariate model`) %>%
  dplyr::select(Metric, `Covariate model`, Predictor,
                `Test statistic`, `p`, `Functional Variance Explained (%)`,
                `Semi-partial r`)

cat(sprintf("Supplementary Table 4: %d rows\n\n", nrow(note4_t4)))
cat("First 14 rows (Strength block):\n")
print(utils::head(note4_t4, 14))

write.csv(note4_t4,
          file.path(tables_dir, "SuppTbl4_per_IV_PFFR_output.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path("tables", "supplement",
                      "SuppTbl4_per_IV_PFFR_output.csv")))


# ----------------------------------------------------------------------------
# 4d  Supplementary Table 5 --- K-basis dimension diagnostics
# ----------------------------------------------------------------------------
# For each of the 6 fully-adjusted PFFR models (Strength, Normalized GE,
# Normalized ACC, SW, Raw GE, Raw ACC), reports the basis dimension (k'),
# effective degrees of freedom (edf), k-index, and p-value for every smooth
# term (Intercept, Group, Age at MRI, eTIV, Sex, Rel_Motion). A non-significant
# p-value combined with k-index >= 1 indicates adequate basis dimension.
#
# Source: refund::pffr.check() printed output, parsed line-by-line. pffr.check
# also draws residual diagnostic plots as a side effect; we route those to a
# scratch PDF device which is discarded.
#
# Yields 36 rows: 6 metrics x 6 smooth terms.
# ----------------------------------------------------------------------------

cat("\n--- 4d: Building Supplementary Table 5 ---\n\n")

# refund needed for pffr.check (mgcv needed transitively for k.check internals)
suppressPackageStartupMessages({
  library(refund)
})

# --- Metric spec: directory + RData prefix, supplement row order ---
n4_t5_metrics <- list(
  list(label = "Strength",
       dir   = "str_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "str"),
  list(label = "Normalized GE",
       dir   = "rand_norm_wei_GE_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "rand_norm_wei_GE"),
  list(label = "Normalized ACC",
       dir   = "rand_norm_wei_ACC_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "rand_norm_wei_ACC"),
  list(label = "Small-Worldness",
       dir   = "rand_norm_wei_SW_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "rand_norm_wei_SW"),
  list(label = "Raw GE",
       dir   = "GE_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "GE"),
  list(label = "Raw ACC",
       dir   = "ACC_FDA_Group_age_at_5y_mri_eTIV_sex_Rel_Motion_11-100_fullsample",
       prefix = "ACC")
)
n4_t5_metric_levels <- vapply(n4_t5_metrics, function(x) x$label, character(1))

# --- Pretty-name map for smooth-term display labels ---
# pffr term names appear as "s(densities.vec):Group", "s(densities.vec)" for
# the intercept, etc. After cleaning the prefix, the raw term names map to:
n4_t5_pretty <- c(
  "Intercept"    = "Intercept",
  "Group"        = "Group (VPT vs FT)",
  "eTIV"         = "Total intracranial volume",
  "age_at_5y_mri" = "Age at MRI",
  "sex"          = "Sex",
  "Rel_Motion"   = "Relative motion"
)
n4_t5_term_order <- c("Intercept", "Group (VPT vs FT)", "Age at MRI",
                      "Total intracranial volume", "Sex", "Relative motion")

# --- Extract k-basis diagnostics from one fully-adjusted RData ---
extract_kcheck <- function(rdata_path, metric_label) {
  if (!file.exists(rdata_path)) {
    warning(sprintf("Missing RData: %s", rdata_path))
    return(NULL)
  }
  e <- new.env()
  load(rdata_path, envir = e)
  if (is.null(e$pffr_fit)) return(NULL)
  
  # The k-index and its p-value are computed by k.check() via an UNSEEDED
  # permutation test, so they vary run-to-run (~+/-0.1) and can flip a
  # threshold-adjacent k-index between a two-decimal print ("1.04") and a bare
  # integer ("1"), shifting column alignment. Seeding makes the diagnostics
  # deterministic and reproducible. k' and edf are deterministic regardless.
  set.seed(5)
  tmp_dev <- tempfile(fileext = ".pdf")
  pdf(file = tmp_dev)
  out_lines <- tryCatch(
    capture.output(suppressWarnings(suppressMessages(
      refund::pffr.check(e$pffr_fit)
    ))),
    error = function(err) {
      message(sprintf("    pffr.check() failed for %s: %s", metric_label, err$message))
      character(0)
    },
    finally = {
      try(dev.off(), silent = TRUE)
    }
  )
  unlink(tmp_dev)
  
  if (length(out_lines) == 0) return(NULL)
  
  # Locate the k-basis table header
  hdr_idx <- grep("^\\s*k'\\s+edf\\s+k-index\\s+p-value", out_lines)
  if (length(hdr_idx) == 0) return(NULL)
  
  body_start <- hdr_idx[1] + 1
  body_lines <- out_lines[body_start:length(out_lines)]
  body_lines <- body_lines[nzchar(trimws(body_lines))]
  
  rows <- lapply(body_lines, function(ln) {
    ln <- trimws(ln)
    # Skip only the legend line ("Signif. codes: ...") and the "---" rule.
    # Do NOT skip on bare "***": data rows with highly significant k-check
    # p-values (e.g. Raw ACC, "<2e-16 ***") carry asterisks and must be kept;
    # the trailing signif code is stripped just below.
    if (grepl("Signif\\. codes", ln) || grepl("^-{3,}$", ln)) return(NULL)
    # Strip any trailing mgcv significance code so it can't be misread as a
    # column or collapse against the p-value.
    ln <- sub("[[:space:]]+(\\*{1,3}|\\.)[[:space:]]*$", "", ln)
    # k.check() prints exactly 4 trailing numeric fields: k' edf k-index p-value.
    # Tolerates bare-integer k-index ("1") and relational p-values ("<2e-16").
    m <- regmatches(
      ln,
      regexec(paste0("^(.*?)\\s+",
                     "([0-9.eE+-]+)\\s+",        # k'
                     "([0-9.eE+-]+)\\s+",        # edf
                     "([0-9.eE+-]+)\\s+",        # k-index
                     "([<>]?[0-9.eE+-]+)\\s*$"), # p-value (may carry < or >)
              ln)
    )[[1]]
    if (length(m) < 6) return(NULL)
    term_name <- trimws(m[2])
    kp        <- suppressWarnings(as.numeric(m[3]))
    edf       <- suppressWarnings(as.numeric(m[4]))
    ki        <- suppressWarnings(as.numeric(m[5]))
    pv        <- trimws(m[6])
    if (!nzchar(term_name) || is.na(kp) || is.na(edf) || is.na(ki)) return(NULL)
    data.frame(
      term       = term_name,
      `k'`       = kp,
      edf        = edf,
      `k-index`  = ki,
      `p-value`  = pv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  
  k_df <- bind_rows(rows)
  k_df$metric <- metric_label
  k_df[, c("metric", "term", "k'", "edf", "k-index", "p-value")]
}

# --- Loop through all 6 fully-adjusted fits ---
t5_blocks <- list()
for (mspec in n4_t5_metrics) {
  rdata_path <- file.path(repo_root, "results", "fda_group_differences",
                          mspec$dir, paste0(mspec$prefix, "_FDA_results.RData"))
  cat(sprintf("  [%s] ", mspec$label))
  k_df <- extract_kcheck(rdata_path, mspec$label)
  if (is.null(k_df)) {
    cat("no output\n"); next
  }
  cat(sprintf("%d terms\n", nrow(k_df)))
  t5_blocks[[length(t5_blocks) + 1]] <- k_df
}

if (length(t5_blocks) == 0) {
  cat("\n  [warning] no k-basis output collected -- check RData paths\n")
} else {
  t5_raw <- bind_rows(t5_blocks)
  
  # Clean term names: "s(densities.vec):Group" -> "Group"; "s(densities.vec)"
  # (no colon) is the functional intercept.
  note4_t5 <- t5_raw %>%
    mutate(
      term_clean = ifelse(grepl(":", term, fixed = TRUE),
                          sub("^[^:]+:", "", term),
                          "Intercept"),
      Term = ifelse(term_clean %in% names(n4_t5_pretty),
                    n4_t5_pretty[term_clean],
                    term_clean),
      `k'`      = as.character(`k'`),     # k' should display as integer (no decimals)
      edf       = sprintf("%.2f", as.numeric(edf)),
      `k-index` = sprintf("%.2f", as.numeric(`k-index`))
    ) %>%
    dplyr::select(Metric = metric, Term, `k'`, edf, `k-index`, `p-value`) %>%
    mutate(
      Metric = factor(Metric, levels = n4_t5_metric_levels),
      Term   = factor(Term,   levels = n4_t5_term_order)
    ) %>%
    arrange(Metric, Term)
  
  cat(sprintf("\nSupplementary Table 5: %d rows\n\n", nrow(note4_t5)))
  cat("First 12 rows (Strength + Normalized GE):\n")
  print(utils::head(note4_t5, 12))
  
  write.csv(note4_t5,
            file.path(tables_dir, "SuppTbl5_kbasis_diagnostics.csv"),
            row.names = FALSE)
  cat(sprintf("\nSaved: %s\n",
              file.path("tables", "supplement",
                        "SuppTbl5_kbasis_diagnostics.csv")))
}


# ============================================================================
# 4e: Supplementary Figures 6-10 - covariate / group beta(density) curves
# ============================================================================
#
# Figs 6-9: one figure per metric (Strength, Normalized GE, Normalized ACC,
#   Small-Worldness), each a 2x2 facet of the four COVARIATE beta(density)
#   curves (TIV, Age at MRI, Sex, Relative motion) from the fully-adjusted
#   PFFR fit, with pointwise 95% bootstrap CI ribbons.
# Fig 10: the GROUP beta(density) curves for the two RAW metrics (Raw GE,
#   Raw ACC), 1x2 facet, rendered in red to distinguish from the normalized
#   metrics' panels.
#
# Source: the same fully-adjusted *_FDA_results.RData used by Sup Tbl 5
#   (n4_t5_metrics), using the n4_t5_pretty label map and the shared color,
#   facet, and theme conventions.
# ============================================================================

cat("\n--- 4e: Building Supplementary Figures 6-10 ---\n\n")

# Covariate panel order (matches figure captions: TIV, age, sex, motion)
covar_order <- c("Total intracranial volume", "Age at MRI", "Sex", "Relative motion")

# --- Capture fitted beta(density) curves from a pffr fit ---------------------
# plot.pffr() returns, invisibly, a list of per-term plot objects each carrying
# $x (density grid) and $fit (fitted beta). Route the plot to a throwaway
# device so nothing is drawn. (Replicates the pipeline's getPlotObject.)
get_plot_object <- function(model) {
  ff <- tempfile(fileext = ".pdf")
  pdf(file = ff)
  po <- plot(model)
  dev.off()
  unlink(ff)
  po
}

# --- Extract beta(d) + bootstrap CI for ALL smooth terms of one fitted model -
# Returns a long-format data.frame: metric, predictor_raw, predictor (pretty),
# density, beta, ci_lower, ci_upper. One block of rows per smooth term.
extract_beta_curves <- function(rdata_path, metric_label) {
  if (!file.exists(rdata_path)) {
    cat(sprintf("  [missing] %s -- skipping\n", basename(rdata_path)))
    return(NULL)
  }
  e <- new.env()
  load(rdata_path, envir = e)
  if (is.null(e$pffr_fit) || is.null(e$bootstrap_coefs)) {
    cat(sprintf("  [no bootstrap] %s\n", basename(rdata_path)))
    return(NULL)
  }
  
  # bootstrap_coefs$smterms names look like "Group(densities)",
  # "eTIV(densities)", etc. Strip the "(...)" suffix for the pretty-map lookup.
  raw_term_names <- names(e$bootstrap_coefs$smterms)
  term_names     <- sub("\\(.*\\)$", "", raw_term_names)
  
  po     <- get_plot_object(e$pffr_fit)
  smList <- e$bootstrap_coefs$smterms
  n_terms <- min(length(po), length(smList), length(term_names))
  
  out_list <- list()
  for (i in seq_len(n_terms)) {
    pred_raw <- term_names[i]
    if (is.na(pred_raw) || pred_raw == "") pred_raw <- "Intercept"
    pred_lbl <- if (pred_raw %in% names(n4_t5_pretty)) n4_t5_pretty[[pred_raw]] else pred_raw
    
    out_list[[i]] <- data.frame(
      metric        = metric_label,
      predictor_raw = pred_raw,
      predictor     = pred_lbl,
      density       = po[[i]]$x,
      beta          = po[[i]]$fit,
      ci_lower      = smList[[i]][["2.5%"]],
      ci_upper      = smList[[i]][["97.5%"]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out_list)
}

# --- Extract beta curves for all 6 fully-adjusted fits -----------------------
beta_all <- list()
for (mspec in n4_t5_metrics) {
  rdata_path <- file.path(repo_root, "results", "fda_group_differences",
                          mspec$dir, paste0(mspec$prefix, "_FDA_results.RData"))
  cat(sprintf("  [%s] ", mspec$label))
  bdf <- extract_beta_curves(rdata_path, mspec$label)
  if (is.null(bdf)) { cat("no output\n"); next }
  beta_all[[mspec$label]] <- bdf
  cat(sprintf("%d terms x %d densities\n",
              length(unique(bdf$predictor)), length(unique(bdf$density))))
}

# --- Figs 6-9: covariate beta(d) curves, one figure per normalized metric ----
make_covariate_beta_fig <- function(metric_label, fig_num) {
  if (!metric_label %in% names(beta_all)) {
    cat(sprintf("  [skip] Fig %d: %s (no beta data)\n", fig_num, metric_label))
    return(invisible(NULL))
  }
  df_one <- beta_all[[metric_label]]
  df_one <- df_one[df_one$predictor %in% covar_order, ]
  if (nrow(df_one) == 0) {
    cat(sprintf("  [skip] Fig %d: %s (no covariate terms)\n", fig_num, metric_label))
    return(invisible(NULL))
  }
  df_one$predictor <- factor(df_one$predictor, levels = covar_order)
  
  p <- ggplot(df_one, aes(x = density, y = beta)) +
    geom_hline(yintercept = 0, color = "grey60",
               linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                fill = "#1F77B4", alpha = 0.25) +
    geom_line(color = "#1F4E79", linewidth = 0.8) +
    facet_wrap(~ predictor, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(20, 100, by = 20),
                       limits = c(density_min, density_max),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Network density (%)",
         y = expression(beta(d)),
         title = sprintf("Supplementary Figure %d. Covariate beta(density) curves: %s",
                         fig_num, metric_label)) +
    theme_supp() +
    theme(strip.text = element_text(size = 10, face = "bold"))
  
  fname <- sprintf("SuppFig%d_covariate_beta_%s",
                   fig_num, gsub("[^A-Za-z0-9]+", "_", metric_label))
  save_fig(p, fname, width_in = 8.5, height_in = 6.5)
}

make_covariate_beta_fig("Strength",        6)
make_covariate_beta_fig("Normalized GE",   7)
make_covariate_beta_fig("Normalized ACC",  8)
make_covariate_beta_fig("Small-Worldness", 9)

# --- Fig 10: Group beta(d) for the two RAW metrics, rendered in red ----------
df_raw <- do.call(rbind, list(beta_all[["Raw GE"]], beta_all[["Raw ACC"]]))
if (!is.null(df_raw) && nrow(df_raw) > 0) {
  df_raw <- df_raw[df_raw$predictor == "Group (VPT vs FT)", ]
  df_raw$metric <- factor(df_raw$metric, levels = c("Raw GE", "Raw ACC"))
}

if (!is.null(df_raw) && nrow(df_raw) > 0) {
  p_fig10 <- ggplot(df_raw, aes(x = density, y = beta)) +
    geom_hline(yintercept = 0, color = "grey60",
               linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                fill = "#D62728", alpha = 0.20) +
    geom_line(color = "#8B1A1A", linewidth = 0.9) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(20, 100, by = 20),
                       limits = c(density_min, density_max),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Network density (%)",
         y = expression(beta(d)),
         title = "Supplementary Figure 10. Group beta(density): Raw GE and Raw ACC") +
    theme_supp() +
    theme(strip.text = element_text(size = 10, face = "bold"))
  save_fig(p_fig10, "SuppFig10_raw_GE_ACC_group_beta",
           width_in = 7.5, height_in = 3.5)
} else {
  cat("  [warning] Fig 10: no raw GE / raw ACC group beta data\n")
}


cat("\n========================================\n")
cat("Note 4 DONE\n")
cat("========================================\n")

# ============================================================================
# ============================================================================
# SUPPLEMENTARY NOTE 5: UNIVARIATE EXPOSURE-METRIC ASSOCIATIONS
#
#   Supplementary Table 6  - Univariate fully-adjusted PFFR results
#                            (7 exposures x 4 primary metrics = 28 rows)
#   Supplementary Figs 11-14 - Univariate beta(density) curves per metric
#                            (7 exposure panels each, 3-3-1 layout)
#
# Each cell is a fully-adjusted univariate PFFR model: one neonatal exposure
# plus the forced covariate set (eTIV, sex, social risk score, age at MRI,
# relative motion), fit within the VPT group only (n=171). All cells use the
# gaussian family, including SW (Methods 2.7.4 tractability caveat).
#
# SuppTbl6 columns:
#   F, p           = the overall-model test for the cell: the anova comparison
#                    of the full model (exposure + forced covariates) vs the
#                    intercept-only model (from summary_results$Overall_*).
#   Significant    = "Y" if the exposure bootstrap beta(d) CI excludes 0 at any
#                    density, else "N".
#   Direction      = sign of mean beta over the significant region.
#   CI from / to   = density range over which the CI excludes 0.
#   % density sig  = fraction of densities where the CI excludes 0.
#
# Source: results/fda_univariate/{metric}/{metric}_univariate_{exposure}_11-100/
#         {metric}_univariate_{exposure}_FDA_results.RData
# Plotting reuses theme_supp(), save_fig(), get_plot_object() (defined in 4e).
# ============================================================================
# ============================================================================

cat("\n========================================\n")
cat("Note 5: Univariate Exposure-Metric Associations\n")
cat("========================================\n\n")

# --- Metric specs (4 primary metrics; outcome folder + RData prefix) ---------
n5_metrics <- list(
  list(label = "Strength",        prefix = "str"),
  list(label = "Normalized GE",   prefix = "rand_norm_wei_GE"),
  list(label = "Normalized ACC",  prefix = "rand_norm_wei_ACC"),
  list(label = "Small-Worldness", prefix = "rand_norm_wei_SW")
)
n5_metric_levels <- vapply(n5_metrics, function(x) x$label, character(1))

# --- Exposure specs (7 candidate neonatal exposures; var name + label) -------
# Order defines panel order in Figs 11-14 and row order within each metric
# block of SuppTbl6 (matches THREE_fda_univariate_runs.R).
n5_exposures <- list(
  list(var = "bpd2",              label = "BPD"),
  list(var = "bw_z",              label = "BWZ"),
  list(var = "ga",                label = "GA"),
  list(var = "globalbrainscore2", label = "GBA"),
  list(var = "anyrop",            label = "ROP"),
  list(var = "sepsis2",           label = "Sepsis"),
  list(var = "dwma_percent",      label = "DWMA")
)
n5_panel_levels <- vapply(n5_exposures, function(x) x$label, character(1))

# --- Locate one univariate RData file in the new fda_univariate layout -------
find_univ_rdata <- function(metric_prefix, exposure_var) {
  run_dir <- sprintf("%s_univariate_%s_11-100", metric_prefix, exposure_var)
  rdata   <- sprintf("%s_univariate_%s_FDA_results.RData", metric_prefix, exposure_var)
  p <- file.path(repo_root, "results", "fda_univariate", metric_prefix, run_dir, rdata)
  if (file.exists(p)) return(p)
  NULL
}

# --- Extract exposure beta(d) curve + SuppTbl6 stats from one cell -----------
# F/p come from the overall-model test (full vs intercept-only) in
# summary_results. Significance/direction/CI/%sig come from the bootstrap beta(d) CI.
extract_univ_beta_and_stats <- function(rdata_path, exposure_var, exposure_label,
                                        metric_label) {
  e <- new.env()
  load(rdata_path, envir = e)
  if (is.null(e$bootstrap_coefs)) return(NULL)
  
  # Find the exposure's smooth term: name pattern "{var}(densities)"
  smterm_names <- names(e$bootstrap_coefs$smterms)
  exp_term <- paste0(exposure_var, "(densities)")
  if (!(exp_term %in% smterm_names)) {
    match_idx <- grep(exposure_var, smterm_names, fixed = TRUE)
    if (length(match_idx) == 0) {
      cat(sprintf("    [warning] no smterm for '%s' in %s\n",
                  exposure_var, basename(rdata_path)))
      return(NULL)
    }
    exp_term <- smterm_names[match_idx[1]]
  }
  
  # beta(d) + bootstrap CI from the plot object (beta) and smterms (CIs),
  # as in extract_beta_curves(). Index by matching the exposure term.
  po       <- get_plot_object(e$pffr_fit)
  raw_names <- sub("\\(.*\\)$", "", smterm_names)
  idx <- which(raw_names == exposure_var)
  if (length(idx) == 0) idx <- grep(exposure_var, smterm_names, fixed = TRUE)
  idx <- idx[1]
  
  sm_df    <- e$bootstrap_coefs$smterms[[exp_term]]
  density  <- po[[idx]]$x
  beta     <- po[[idx]]$fit
  ci_lower <- sm_df[["2.5%"]]
  ci_upper <- sm_df[["97.5%"]]
  
  # Significance pattern from the bootstrap CI (excludes zero where lower>0 or upper<0)
  sig_mask <- (ci_lower > 0 & ci_upper > 0) | (ci_lower < 0 & ci_upper < 0)
  any_sig  <- any(sig_mask, na.rm = TRUE)
  if (any_sig) {
    sig_d     <- density[sig_mask]
    sig_b     <- beta[sig_mask]
    direction <- if (mean(sig_b, na.rm = TRUE) > 0) "Positive" else "Negative"
    sig_min   <- min(sig_d, na.rm = TRUE)
    sig_max   <- max(sig_d, na.rm = TRUE)
    pct_sig   <- 100 * sum(sig_mask, na.rm = TRUE) / sum(!is.na(sig_mask))
    sig_str   <- "Y"
  } else {
    direction <- NA_character_; sig_min <- NA_real_; sig_max <- NA_real_
    pct_sig <- 0; sig_str <- "N"
  }
  
  # F / p: overall-model test for this cell -- the anova comparison of the full
  # model (exposure + forced covariates) against the intercept-only model,
  # computed in the pipeline and stored in summary_results. This is the overall
  # model F (gaussian cells), identical in definition to the group-difference
  # overall test.
  f_stat <- p_val <- NA_real_
  if (!is.null(e$summary_results)) {
    sr <- e$summary_results
    if ("Overall_Test_Statistic" %in% names(sr)) f_stat <- sr$Overall_Test_Statistic[1]
    if ("Overall_p"              %in% names(sr)) p_val  <- sr$Overall_p[1]
  }
  
  curve_df <- data.frame(
    metric = metric_label, exposure = exposure_label,
    density = density, beta = beta,
    ci_lower = ci_lower, ci_upper = ci_upper,
    stringsAsFactors = FALSE
  )
  stats_row <- data.frame(
    Metric = metric_label, Exposure = exposure_label,
    `F` = f_stat, `p` = p_val,
    Significant = sig_str, Direction = direction,
    `CI excludes 0 from (%)` = sig_min, `to (%)` = sig_max,
    `% density significant` = pct_sig,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  list(curve = curve_df, stats = stats_row)
}

# --- Loop over all 28 (metric x exposure) cells ------------------------------
cat("Extracting univariate fully-adjusted results...\n\n")
curves_all <- list()
stats_all  <- list()
for (mspec in n5_metrics) {
  cat(sprintf("[%s]\n", mspec$label))
  for (espec in n5_exposures) {
    rdata_path <- find_univ_rdata(mspec$prefix, espec$var)
    if (is.null(rdata_path)) {
      cat(sprintf("  [missing] %s x %s\n", mspec$label, espec$label)); next
    }
    res <- extract_univ_beta_and_stats(rdata_path, espec$var, espec$label, mspec$label)
    if (is.null(res)) next
    key <- paste(mspec$label, espec$label, sep = "::")
    curves_all[[key]] <- res$curve
    stats_all[[key]]  <- res$stats
    sig_chr <- res$stats$Significant
    cat(sprintf("  %-7s significant=%s%s\n", espec$label, sig_chr,
                if (sig_chr == "Y")
                  sprintf("  (%s, %.1f-%.1f%%)", res$stats$Direction,
                          res$stats$`CI excludes 0 from (%)`, res$stats$`to (%)`)
                else ""))
  }
  cat("\n")
}

# --- Supplementary Table 6 ---------------------------------------------------
cat("--- Building Supplementary Table 6 ---\n\n")
if (length(stats_all) > 0) {
  fmt_p6 <- function(p) ifelse(is.na(p), NA_character_,
                               ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  n5_table6 <- do.call(rbind, stats_all)
  n5_table6$Metric   <- factor(n5_table6$Metric,   levels = n5_metric_levels)
  n5_table6$Exposure <- factor(n5_table6$Exposure, levels = n5_panel_levels)
  n5_table6 <- n5_table6[order(n5_table6$Metric, n5_table6$Exposure), ]
  n5_table6$`F` <- ifelse(is.na(n5_table6$`F`), NA_character_, sprintf("%.2f", n5_table6$`F`))
  n5_table6$`p` <- fmt_p6(n5_table6$`p`)
  n5_table6$`CI excludes 0 from (%)` <- ifelse(is.na(n5_table6$`CI excludes 0 from (%)`),
                                               NA_character_, sprintf("%.2f", n5_table6$`CI excludes 0 from (%)`))
  n5_table6$`to (%)` <- ifelse(is.na(n5_table6$`to (%)`),
                               NA_character_, sprintf("%.2f", n5_table6$`to (%)`))
  n5_table6$`% density significant` <- sprintf("%.1f", as.numeric(n5_table6$`% density significant`))
  
  cat(sprintf("Supplementary Table 6: %d rows\n\n", nrow(n5_table6)))
  print(utils::head(n5_table6, 14), row.names = FALSE)
  
  write.csv(n5_table6,
            file.path(tables_dir, "SuppTbl6_univariate_exposure_metric_results.csv"),
            row.names = FALSE)
  cat(sprintf("\nSaved: %s\n",
              file.path("tables", "supplement",
                        "SuppTbl6_univariate_exposure_metric_results.csv")))
} else {
  cat("  [warning] no univariate stats collected\n")
}

# --- Supplementary Figures 11-14: per-metric univariate beta(d) panels -------
# Layout: 3 - 3 - 1 (BPD/BWZ/GA / GBA/ROP/Sepsis / DWMA centered).
cat("\n--- Building Supplementary Figures 11-14 ---\n\n")

make_univ_metric_fig <- function(metric_label, fig_num) {
  metric_curves <- do.call(rbind,
                           Filter(function(df) df$metric[1] == metric_label, curves_all))
  if (is.null(metric_curves) || nrow(metric_curves) == 0) {
    cat(sprintf("  [skip] Fig %d: %s (no curves)\n", fig_num, metric_label))
    return(invisible(NULL))
  }
  metric_curves$exposure <- factor(metric_curves$exposure, levels = n5_panel_levels)
  
  build_panel <- function(exp_lbl) {
    df <- subset(metric_curves, exposure == exp_lbl)
    if (nrow(df) == 0) return(patchwork::plot_spacer())
    ggplot(df, aes(x = density, y = beta)) +
      geom_hline(yintercept = 0, color = "grey60",
                 linetype = "dashed", linewidth = 0.4) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                  fill = "#1F77B4", alpha = 0.25) +
      geom_line(color = "#1F4E79", linewidth = 0.8) +
      scale_x_continuous(breaks = seq(20, 100, by = 20),
                         limits = c(density_min, density_max),
                         expand = expansion(mult = c(0.01, 0.01))) +
      labs(x = NULL, y = NULL, title = exp_lbl) +
      theme_supp() +
      theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 9))
  }
  
  panels <- lapply(n5_panel_levels, build_panel)
  names(panels) <- n5_panel_levels
  row1 <- panels[["BPD"]] | panels[["BWZ"]] | panels[["GA"]]
  row2 <- panels[["GBA"]] | panels[["ROP"]] | panels[["Sepsis"]]
  row3 <- patchwork::plot_spacer() | panels[["DWMA"]] | patchwork::plot_spacer()
  combined <- (row1 / row2 / row3) +
    patchwork::plot_annotation(
      title   = sprintf("Supplementary Figure %d. Univariate beta(density) curves: %s",
                        fig_num, metric_label),
      caption = expression("y-axis = " ~ beta(d) * "; x-axis = network density (%)")
    ) &
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  fname <- sprintf("SuppFig%d_univariate_beta_%s",
                   fig_num, gsub("[^A-Za-z0-9]+", "_", metric_label))
  save_fig(combined, fname, width_in = 9.0, height_in = 8.5)
}

make_univ_metric_fig("Strength",        11)
make_univ_metric_fig("Normalized GE",   12)
make_univ_metric_fig("Normalized ACC",  13)
make_univ_metric_fig("Small-Worldness", 14)

cat("\n========================================\n")
cat("Note 5: DONE\n")
cat("========================================\n")


# ============================================================================
# ============================================================================
# SUPPLEMENTARY NOTE 6: ACC STABILITY SELECTION SENSITIVITY ANALYSES
#
#   Supplementary Figure 15 - ACC stability selection, 4 highest-GBA VPT
#                             participants removed (high_gba_rem branch).
#                             GBA remains the only stable exposure (~88%).
#   Supplementary Figure 16 - ACC stability selection, GBA recoded binary
#                             (globalcatmod >=8 vs <8; gba_binary branch).
#                             GBA (binary) remains the only stable exposure (~91%).
#
# Each figure has two panels, identical in structure to main Figure 3:
#   A) selection-frequency bars across the candidate exposures, 70% threshold line
#   B) GBA beta(density) curve from the post-selection PFFR fit, with bootstrap CI
#
# Source (from FOUR_fda_stability_selection_runs.R):
#   results/fda_stability_selection/ACC_sensitivity/
#       rand_norm_wei_ACC_stabsel_11-100_high_gba_rem/rand_norm_wei_ACC_stabsel_results.RData
#       rand_norm_wei_ACC_stabsel_11-100_gba_binary/  rand_norm_wei_ACC_stabsel_results.RData
#
# In the gba_binary branch the GBA variable is "globalcatmod" (not
# "globalbrainscore2"); exposure_labels maps both to "GBA".
#
# These panels follow the same conventions as main Figure 3.
# ============================================================================
# ============================================================================

cat("\n========================================\n")
cat("Note 6: ACC Stability Selection Sensitivity\n")
cat("========================================\n\n")

# --- Stability-selection helpers ---------------------------------------------
stability_threshold <- 0.70

# Exposure display labels. Includes globalcatmod (binary GBA branch) -> "GBA".
n6_exposure_labels <- c(
  bpd2              = "BPD",
  bw_z              = "BWZ",
  ga                = "GA",
  globalbrainscore2 = "GBA",
  globalcatmod      = "GBA",
  anyrop            = "ROP",
  sepsis2           = "Sepsis"
)
label_exposure6 <- function(raw_name) {
  if (raw_name %in% names(n6_exposure_labels)) n6_exposure_labels[[raw_name]] else raw_name
}

# Find a smooth term by predictor name in bootstrap_coefs$smterms.
find_term_index6 <- function(term_names, predictor) {
  bare <- sub("\\(.*\\)$", "", term_names)
  idx <- which(bare == predictor)
  if (length(idx) == 0) idx <- grep(paste0("^", predictor), term_names)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

# CIs list from bootstrap_coefs$smterms.
get_cis_list6 <- function(coefboot_bs) {
  smList <- coefboot_bs$smterms
  out <- vector("list", length(smList))
  for (i in seq_along(smList)) {
    out[[i]] <- list(ci_lower = smList[[i]][["2.5%"]],
                     ci_upper = smList[[i]][["97.5%"]])
  }
  out
}

# Extract beta(d) + bootstrap CI for one predictor from a fitted model in env.
# model_obj_name = "pffr_final" for the post-selection stab-sel fit.
extract_beta_curve6 <- function(env, predictor, model_obj_name = "pffr_final") {
  model <- env[[model_obj_name]]
  if (is.null(model)) stop(sprintf("No %s in environment", model_obj_name))
  if (is.null(env$bootstrap_coefs)) stop("No bootstrap_coefs in environment")
  term_names <- names(env$bootstrap_coefs$smterms)
  idx <- find_term_index6(term_names, predictor)
  if (is.na(idx)) stop(sprintf("Predictor '%s' not found among: %s",
                               predictor, paste(term_names, collapse = ", ")))
  po  <- get_plot_object(model)   # defined in Note 4e
  CIs <- get_cis_list6(env$bootstrap_coefs)
  data.frame(density = po[[idx]]$x, beta = po[[idx]]$fit,
             ci_lower = CIs[[idx]]$ci_lower, ci_upper = CIs[[idx]]$ci_upper)
}

# Effect direction from the bootstrap CI (sign where CI excludes 0; else peak).
exposure_direction6 <- function(beta_df) {
  if (nrow(beta_df) == 0) return(NA_character_)
  excludes <- (beta_df$ci_lower > 0) | (beta_df$ci_upper < 0)
  if (any(excludes, na.rm = TRUE)) {
    sign_val <- sign(mean(beta_df$beta[excludes], na.rm = TRUE))
  } else {
    sign_val <- sign(beta_df$beta[which.max(abs(beta_df$beta))])
  }
  if (is.na(sign_val) || sign_val == 0) return(NA_character_)
  if (sign_val < 0) "Negative" else "Positive"
}

# Density range over which the CI excludes zero (for console reporting).
ci_exclusion_range6 <- function(beta_df) {
  if (nrow(beta_df) == 0) return(NA_character_)
  excludes <- (beta_df$ci_lower > 0) | (beta_df$ci_upper < 0)
  if (!any(excludes, na.rm = TRUE)) return("ns")
  d <- beta_df$density
  sprintf("%.1f-%.1f%%", min(d[excludes], na.rm = TRUE), max(d[excludes], na.rm = TRUE))
}

# Panel A: selection-frequency bars.
build_stabsel_freq_panel <- function(selection_freq, panel_letter, panel_title) {
  if (is.null(selection_freq) || length(selection_freq) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("%s. %s\n(PENDING)", panel_letter, panel_title)) +
             theme_void())
  }
  freq_df <- data.frame(
    exposure_display = vapply(names(selection_freq), label_exposure6, ""),
    frequency_pct    = as.numeric(selection_freq) * 100,
    stringsAsFactors = FALSE
  )
  freq_df <- freq_df[order(freq_df$frequency_pct, decreasing = TRUE), ]
  freq_df$exposure_display <- factor(freq_df$exposure_display,
                                     levels = freq_df$exposure_display)
  freq_df$above <- freq_df$frequency_pct >= stability_threshold * 100
  freq_df$label <- sprintf("%d%%", as.integer(round(freq_df$frequency_pct)))
  ggplot(freq_df, aes(x = exposure_display, y = frequency_pct, fill = above)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = stability_threshold * 100,
               lty = "dashed", color = "grey50", linewidth = 0.6) +
    geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("FALSE" = "#5b89c7", "TRUE" = "#1a237e"), guide = "none") +
    scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
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

# Panel B: GBA beta(d) curve.
build_stabsel_beta_panel <- function(beta_df, panel_letter, panel_title) {
  if (is.null(beta_df) || nrow(beta_df) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("%s. %s\n(PENDING)", panel_letter, panel_title)) +
             theme_void())
  }
  ggplot(beta_df, aes(x = density, y = beta)) +
    geom_hline(yintercept = 0, lty = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") +
    geom_line(linewidth = 1, color = "darkblue") +
    ylab(expression(beta(d))) + xlab("Density (%)") +
    ggtitle(sprintf("%s. %s", panel_letter, panel_title)) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}

# --- Sensitivity-branch specs ------------------------------------------------
# Each branch: the ACC_sensitivity subdir, the GBA term name in that branch,
# the figure number, and a title tag for panel A.
n6_branches <- list(
  list(branch    = "high_gba_rem",
       subdir    = "rand_norm_wei_ACC_stabsel_11-100_high_gba_rem",
       gba_term  = "globalbrainscore2",
       fig_num   = 15,
       freq_tag  = "ACC stability selection (high-GBA removed)",
       beta_tag  = "GBA effect on ACC (high-GBA removed)"),
  list(branch    = "gba_binary",
       subdir    = "rand_norm_wei_ACC_stabsel_11-100_gba_binary",
       gba_term  = "globalcatmod",
       fig_num   = 16,
       freq_tag  = "ACC stability selection (GBA binary)",
       beta_tag  = "GBA (binary) effect on ACC")
)

stabsel_root6 <- file.path(repo_root, "results", "fda_stability_selection", "ACC_sensitivity")
rdata_name6   <- "rand_norm_wei_ACC_stabsel_results.RData"

for (b in n6_branches) {
  cat(sprintf("--- Supplementary Figure %d (%s) ---\n", b$fig_num, b$branch))
  rdata_path <- file.path(stabsel_root6, b$subdir, rdata_name6)
  if (!file.exists(rdata_path)) {
    cat(sprintf("  [missing] %s -- skipping Fig %d\n", rdata_path, b$fig_num)); next
  }
  env <- new.env(); load(rdata_path, envir = env)
  
  selfreq <- env$selection_freq
  if (is.null(selfreq)) { cat("  [skip] no selection_freq present\n"); next }
  
  # Report top exposure / stable set
  ord  <- order(selfreq, decreasing = TRUE)
  top1 <- names(selfreq)[ord][1]
  cat(sprintf("  Top exposure: %s (%d%%)\n",
              label_exposure6(top1), as.integer(round(selfreq[ord][1] * 100))))
  stable <- env$stable_exposures; if (is.null(stable)) stable <- character(0)
  cat(sprintf("  Stable exposure(s): %s\n",
              if (length(stable)) paste(vapply(stable, label_exposure6, ""), collapse = ", ") else "none"))
  
  # Panel A
  freq_panel <- build_stabsel_freq_panel(selfreq, "A", b$freq_tag)
  
  # Panel B: GBA beta(d) from the post-selection fit
  beta_panel <- NULL
  if (!is.null(env$bootstrap_coefs)) {
    beta_df <- tryCatch(
      extract_beta_curve6(env, b$gba_term, model_obj_name = "pffr_final"),
      error = function(e) { cat(sprintf("  [warn] beta extract failed: %s\n", conditionMessage(e))); NULL })
    if (!is.null(beta_df)) {
      cat(sprintf("  GBA beta(d): direction %s, CI %s\n",
                  exposure_direction6(beta_df), ci_exclusion_range6(beta_df)))
      beta_panel <- build_stabsel_beta_panel(beta_df, "B", b$beta_tag)
    }
  } else {
    cat("  [warn] bootstrap_coefs is NULL; panel B omitted\n")
  }
  
  fig <- if (!is.null(beta_panel)) (freq_panel | beta_panel) else freq_panel
  save_fig(fig, sprintf("SuppFig%d_ACC_stabsel_%s", b$fig_num, b$branch),
           width_in = 10, height_in = 4.5)
  cat("\n")
}

cat("\n========================================\n")
cat("Note 6: DONE\n")
cat("========================================\n")


# ============================================================================
# ============================================================================
# SUPPLEMENTARY NOTE 7: SHARED-SEVERITY ASSESSMENT FRAMEWORK
#
#   Supplementary Figure 17 - Inter-exposure correlation heatmap
#   Supplementary Figure 18 - PCA scree (A) + PC1/PC2 loadings (B)
#   Supplementary Figure 19 - 8-panel scatterplots: exposure PC1 vs FPC1/FPC2
#                             per metric
#   Supplementary Figure 20 - 3-panel permutation results
#                             (A=global, B=incremental heatmap, C=latent severity)
#   Supplementary Table 7   - Permutation test results (3 tests x metrics)
#
# Reads committed outputs from the standalone analysis script
# (note7_shared_severity_analysis.R), in results/post_lasso_shared_severity/.
# This builder only styles; it performs no analysis. Fig 19 reads the
# precomputed per-subject PC1_FPC_scores.csv (no in-builder FPCA).
#
# Reads the precomputed per-subject PC1/FPC scores to render the panels.
# ============================================================================
# ============================================================================

cat("\n========================================\n")
cat("Note 7: Shared-Severity Assessment Framework\n")
cat("========================================\n\n")

n7_root <- file.path(repo_root, "results", "post_lasso_shared_severity")

# Map the analysis script's display labels -> manuscript-consistent labels.
n7_relabel <- c(
  "BPD"            = "BPD",
  "Birth Weight Z" = "BWZ",
  "Gestational Age" = "GA",
  "GBS2"           = "GBA",
  "ROP"            = "ROP",
  "Sepsis"         = "Sepsis"
)
n7_exposure_levels <- c("BPD", "BWZ", "GA", "GBA", "ROP", "Sepsis")
n7_metric_levels   <- c("Strength", "Normalized GE", "Normalized ACC", "Small-Worldness")

# ----------------------------------------------------------------------------
# Supplementary Figure 17: Inter-exposure correlation heatmap
# ----------------------------------------------------------------------------
cat("--- Building Supplementary Figure 17 (correlation heatmap) ---\n")

corr_csv <- file.path(n7_root, "exposure_correlation_matrix.csv")
if (file.exists(corr_csv)) {
  corr_df_raw <- read.csv(corr_csv, row.names = 1, check.names = FALSE)
  rn <- rownames(corr_df_raw); cn <- colnames(corr_df_raw)
  rownames(corr_df_raw) <- ifelse(rn %in% names(n7_relabel), n7_relabel[rn], rn)
  colnames(corr_df_raw) <- ifelse(cn %in% names(n7_relabel), n7_relabel[cn], cn)
  if (all(n7_exposure_levels %in% rownames(corr_df_raw))) {
    corr_df_raw <- corr_df_raw[n7_exposure_levels, n7_exposure_levels]
  }
  corr_long <- corr_df_raw %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "r") %>%
    mutate(Var1 = factor(Var1, levels = rownames(corr_df_raw)),
           Var2 = factor(Var2, levels = colnames(corr_df_raw)))
  p_fig17 <- ggplot(corr_long, aes(x = Var2, y = Var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3.2, color = "black") +
    scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#D62728",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = "Supplementary Figure 17. Inter-exposure correlation matrix") +
    theme_supp() +
    theme(axis.text.x = element_text(angle = 0))
  save_fig(p_fig17, "SuppFig17_exposure_correlation_heatmap",
           width_in = 6.5, height_in = 5.5)
} else {
  cat(sprintf("  [missing] %s\n", corr_csv))
}

# ----------------------------------------------------------------------------
# Supplementary Figure 18: PCA scree (A) + PC1/PC2 loadings (B)
# ----------------------------------------------------------------------------
cat("\n--- Building Supplementary Figure 18 (PCA scree + loadings) ---\n")

pca_var_csv  <- file.path(n7_root, "exposure_pca_variance.csv")
pca_load_csv <- file.path(n7_root, "exposure_pca_loadings.csv")
if (file.exists(pca_var_csv) && file.exists(pca_load_csv)) {
  pca_var  <- read.csv(pca_var_csv, stringsAsFactors = FALSE)
  pca_load <- read.csv(pca_load_csv, stringsAsFactors = FALSE)
  pca_load$Exposure <- ifelse(pca_load$Exposure %in% names(n7_relabel),
                              n7_relabel[pca_load$Exposure], pca_load$Exposure)
  
  # Panel A: scree (bars = per-PC variance; dashed red line = cumulative)
  p18_A <- ggplot(pca_var, aes(x = PC, y = Variance_Explained, group = 1)) +
    geom_col(fill = "#1F77B4", alpha = 0.7, width = 0.6) +
    geom_line(aes(y = Cumulative), color = "#D62728", linewidth = 0.8,
              linetype = "dashed") +
    geom_point(aes(y = Cumulative), color = "#D62728", size = 2) +
    geom_text(aes(label = sprintf("%.0f%%", Variance_Explained)),
              vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(x = "Principal component", y = "Variance explained (%)",
         title = "A. Variance explained per PC",
         subtitle = "Bars = per-PC; dashed line = cumulative") +
    theme_supp() +
    theme(plot.title = element_text(face = "bold", hjust = 0))
  
  # Panel B: PC1 + PC2 loadings
  if (all(c("PC1", "PC2") %in% names(pca_load))) {
    load_long <- pca_load %>%
      dplyr::select(Exposure, PC1, PC2) %>%
      pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
      mutate(Exposure = factor(Exposure, levels = n7_exposure_levels))
    p18_B <- ggplot(load_long, aes(x = Exposure, y = loading, fill = PC)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
      scale_fill_manual(values = c("PC1" = "#1F77B4", "PC2" = "#FF7F0E")) +
      labs(x = NULL, y = "Loading", title = "B. PC1 and PC2 loadings",
           subtitle = "PC1 = latent severity axis", fill = NULL) +
      theme_supp() +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", hjust = 0))
  } else {
    p18_B <- patchwork::plot_spacer()
  }
  
  combined_18 <- (p18_A | p18_B) +
    patchwork::plot_annotation(
      title = "Supplementary Figure 18. Principal component analysis of exposure set",
      theme = theme(plot.title = element_text(size = 12, face = "bold")))
  save_fig(combined_18, "SuppFig18_exposure_PCA", width_in = 9.5, height_in = 4.5)
} else {
  cat("  [missing] PCA CSVs\n")
}

# ----------------------------------------------------------------------------
# Supplementary Figure 19: 8-panel PC1 vs FPC1/FPC2 scatterplots per metric
# Reads precomputed PC1_FPC_scores.csv (no in-builder FPCA).
# ----------------------------------------------------------------------------
cat("\n--- Building Supplementary Figure 19 (PC1 vs FPC scatterplots) ---\n")

scores_csv <- file.path(n7_root, "PC1_FPC_scores.csv")
if (file.exists(scores_csv)) {
  # The analysis script writes metric_label in its own convention
  # ("Clustering (ACC)", etc.); map to manuscript labels before factoring,
  # exactly as relabel_metric7() does for the permutation data. Without this,
  # factor(levels = n7_metric_levels) coerces every value to NA and the
  # per-metric panel filter matches zero rows.
  n7_scores_metric_relabel <- c(
    "Clustering (ACC)"       = "Normalized ACC",
    "Global Efficiency (GE)" = "Normalized GE",
    "Strength"               = "Strength",
    "Small-Worldness (SW)"   = "Small-Worldness"
  )
  n7_scores_long <- read.csv(scores_csv, stringsAsFactors = FALSE) %>%
    mutate(
      metric_label = ifelse(metric_label %in% names(n7_scores_metric_relabel),
                            n7_scores_metric_relabel[metric_label], metric_label),
      metric_label = factor(metric_label, levels = n7_metric_levels)
    )
  
  build_scatter_panel <- function(df, fpc_col, metric_lab, fpc_letter) {
    df_p <- df %>% filter(metric_label == metric_lab,
                          !is.na(.data[[fpc_col]]), !is.na(PC1))
    # Guard: cor.test needs >= 3 finite pairs. If a panel has too few rows
    # (e.g. a metric absent from the scores CSV), render an empty placeholder
    # rather than erroring the whole figure.
    if (nrow(df_p) < 3) {
      cat(sprintf("  [warn] %s / %s: %d finite pairs (<3); placeholder panel\n",
                  metric_lab, fpc_col, nrow(df_p)))
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = sprintf("%s. %s vs. PC1\n%s\n(insufficient data)",
                                        fpc_letter, fpc_col, metric_lab)) +
               theme_void())
    }
    ct    <- suppressWarnings(cor.test(df_p$PC1, df_p[[fpc_col]], method = "pearson"))
    r_val <- ct$estimate; p_val <- ct$p.value; n_n <- nrow(df_p)
    p_lab <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
    r_lab <- sprintf("r = %.2f, %s (n = %d)", r_val, p_lab, n_n)
    ggplot(df_p, aes(x = PC1, y = .data[[fpc_col]])) +
      geom_point(alpha = 0.55, color = "#1F77B4", size = 1.6) +
      geom_smooth(method = "lm", se = TRUE, color = "#D62728",
                  fill = "#D62728", alpha = 0.18, linewidth = 0.7, formula = y ~ x) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
               label = r_lab, size = 3, color = "grey20") +
      labs(x = "Exposure PC1 (latent severity)", y = fpc_col,
           title = sprintf("%s. %s vs. PC1", fpc_letter, fpc_col),
           subtitle = metric_lab) +
      theme_supp() +
      theme(plot.title = element_text(face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 9, color = "grey30"))
  }
  
  panel_letters <- c("A", "B", "C", "D", "E", "F", "G", "H")
  panels <- list(); li <- 1
  for (ml in n7_metric_levels) {
    panels[[length(panels) + 1]] <- build_scatter_panel(n7_scores_long, "FPC1", ml, panel_letters[li]); li <- li + 1
    panels[[length(panels) + 1]] <- build_scatter_panel(n7_scores_long, "FPC2", ml, panel_letters[li]); li <- li + 1
  }
  combined_19 <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "Supplementary Figure 19. Exposure PC1 vs FPC scores per graph metric",
      theme = theme(plot.title = element_text(size = 12, face = "bold")))
  save_fig(combined_19, "SuppFig19_PC1_vs_FPC_scatterplots", width_in = 9.5, height_in = 13)
} else {
  cat(sprintf("  [missing] %s\n", scores_csv))
}

# ----------------------------------------------------------------------------
# Supplementary Figure 20: 3-panel permutation results
# A = global, B = incremental heatmap, C = latent severity
# ----------------------------------------------------------------------------
cat("\n--- Building Supplementary Figure 20 (permutation results 3-panel) ---\n")

global_csv <- file.path(n7_root, "perm_global_exposure_set_by_metric.csv")
inc_csv    <- file.path(n7_root, "perm_incremental_exposure_by_metric.csv")
latent_csv <- file.path(n7_root, "perm_latent_severity_by_metric.csv")
global_df <- if (file.exists(global_csv)) read.csv(global_csv, stringsAsFactors = FALSE) else NULL
inc_df    <- if (file.exists(inc_csv))    read.csv(inc_csv,    stringsAsFactors = FALSE) else NULL
latent_df <- if (file.exists(latent_csv)) read.csv(latent_csv, stringsAsFactors = FALSE) else NULL

relabel_metric7 <- function(df) {
  if (is.null(df) || !("metric_label" %in% names(df))) return(df)
  mr <- c("Clustering (ACC)" = "Normalized ACC",
          "Global Efficiency (GE)" = "Normalized GE",
          "Strength" = "Strength",
          "Small-Worldness (SW)" = "Small-Worldness")
  df$metric_label <- ifelse(df$metric_label %in% names(mr), mr[df$metric_label], df$metric_label)
  df$metric_label <- factor(df$metric_label, levels = n7_metric_levels)
  df
}
global_df <- relabel_metric7(global_df)
inc_df    <- relabel_metric7(inc_df)
latent_df <- relabel_metric7(latent_df)
if (!is.null(inc_df) && "exposure_label" %in% names(inc_df)) {
  inc_df$exposure_label <- ifelse(inc_df$exposure_label %in% names(n7_relabel),
                                  n7_relabel[inc_df$exposure_label], inc_df$exposure_label)
  inc_df$exposure_label <- factor(inc_df$exposure_label, levels = n7_exposure_levels)
}

build_perm_panel <- function(df, value_col, panel_letter, panel_title, panel_subtitle = NULL) {
  if (is.null(df) || !(value_col %in% names(df))) return(patchwork::plot_spacer())
  df$p_val <- pmax(as.numeric(df[[value_col]]), 1e-4)
  df$logp  <- -log10(df$p_val)
  df$p_label <- ifelse(df$p_val < 0.001, "<0.001", sprintf("%.3f", df$p_val))
  ggplot(df, aes(x = reorder(metric_label, logp), y = logp)) +
    geom_col(fill = "#1F77B4", alpha = 0.85, width = 0.6) +
    geom_text(aes(label = p_label), hjust = -0.15, size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "#D62728", linewidth = 0.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = NULL, y = expression(-log[10](p[perm])),
         title = sprintf("%s. %s", panel_letter, panel_title), subtitle = panel_subtitle) +
    theme_supp() +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

p20_A <- build_perm_panel(global_df, "p_global_perm", "A", "Global exposure set",
                          "Joint association of all exposures with FPCA scores")
if (!is.null(inc_df) && "p_inc_perm" %in% names(inc_df)) {
  inc_plot_df <- inc_df %>%
    mutate(p_val = pmax(as.numeric(p_inc_perm), 1e-4),
           logp  = -log10(p_val),
           p_label = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)))
  p20_B <- ggplot(inc_plot_df, aes(x = exposure_label, y = metric_label, fill = logp)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = p_label), size = 3) +
    scale_fill_gradient(low = "white", high = "#1F4E79", name = expression(-log[10](p))) +
    labs(x = NULL, y = NULL, title = "B. Incremental (per-exposure)",
         subtitle = "Each exposure's unique contribution beyond the others") +
    theme_supp() +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          axis.text.x = element_text(angle = 0))
} else {
  p20_B <- patchwork::plot_spacer()
}
p20_C <- build_perm_panel(latent_df, "p_severity_perm", "C", "Latent severity (PC1)",
                          "Single-axis exposure summary")

combined_20perm <- (p20_A | p20_B) / (p20_C | patchwork::plot_spacer()) +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    title = "Supplementary Figure 20. Permutation-based exposure-metric association tests",
    theme = theme(plot.title = element_text(size = 12, face = "bold")))
save_fig(combined_20perm, "SuppFig20_permutation_results", width_in = 11, height_in = 8)

# ----------------------------------------------------------------------------
# Supplementary Table 7: Combined permutation results
# ----------------------------------------------------------------------------
cat("\n--- Building Supplementary Table 7 (permutation results combined) ---\n")

fmt_perm_p <- function(p) ifelse(is.na(p), NA_character_,
                                 ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
n7_table7_rows <- list()
if (!is.null(global_df)) {
  n7_table7_rows[["global"]] <- global_df %>%
    transmute(Test = "Global exposure set", Metric = metric_label, Exposure = "All",
              `K (FPCs)` = if ("K" %in% names(global_df)) K else NA_integer_,
              `cumPVE (%)` = if ("achieved_cumPVE" %in% names(global_df))
                sprintf("%.1f", as.numeric(achieved_cumPVE) * 100) else NA_character_,
              `p (perm)` = fmt_perm_p(as.numeric(p_global_perm)))
}
if (!is.null(inc_df)) {
  n7_table7_rows[["incremental"]] <- inc_df %>%
    transmute(Test = "Incremental (per-exposure)", Metric = metric_label, Exposure = exposure_label,
              `K (FPCs)` = if ("K" %in% names(inc_df)) K else NA_integer_,
              `cumPVE (%)` = if ("achieved_cumPVE" %in% names(inc_df))
                sprintf("%.1f", as.numeric(achieved_cumPVE) * 100) else NA_character_,
              `p (perm)` = fmt_perm_p(as.numeric(p_inc_perm)))
}
if (!is.null(latent_df)) {
  n7_table7_rows[["latent"]] <- latent_df %>%
    transmute(Test = "Latent severity (PC1)", Metric = metric_label, Exposure = "PC1",
              `K (FPCs)` = if ("K" %in% names(latent_df)) K else NA_integer_,
              `cumPVE (%)` = if ("achieved_cumPVE" %in% names(latent_df))
                sprintf("%.1f", as.numeric(achieved_cumPVE) * 100) else NA_character_,
              `p (perm)` = fmt_perm_p(as.numeric(p_severity_perm)))
}
if (length(n7_table7_rows) > 0) {
  n7_table7 <- bind_rows(n7_table7_rows) %>%
    mutate(Test = factor(Test, levels = c("Global exposure set",
                                          "Incremental (per-exposure)",
                                          "Latent severity (PC1)"))) %>%
    arrange(Test, Metric, Exposure)
  cat(sprintf("Supplementary Table 7: %d rows\n\n", nrow(n7_table7)))
  print(utils::head(n7_table7, 12), row.names = FALSE)
  write.csv(n7_table7,
            file.path(tables_dir, "SuppTbl7_permutation_results.csv"), row.names = FALSE)
  cat(sprintf("\nSaved: %s\n",
              file.path("tables", "supplement", "SuppTbl7_permutation_results.csv")))
} else {
  cat("  [warning] no permutation results collected\n")
}

cat("\n========================================\n")
cat("Note 7: DONE\n")
cat("========================================\n")