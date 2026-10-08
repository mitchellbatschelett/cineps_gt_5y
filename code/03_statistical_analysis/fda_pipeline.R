################################################################################
#
#   fda_pipeline.R
#
#   Function-on-scalar regression (pffr) of one graph-theory metric(density)
#   on a set of scalar predictors. One call = one model. Shared by the
#   group-difference and univariate analyses: driven by
#   TWO_fda_group_difference_runs.R or THREE_fda_univariate_runs.R, which set
#   the BATCH_* variables below and source this file.
#
#   Model specification:
#     - pffr, bs.yindex = ps, k = 20, m = c(2, 1), density 11-100%
#     - continuous predictors z-scored; binary predictors not scaled
#     - overall LR test (varying vs constant coefficients) + per-predictor LR
#     - F test for gaussian, Chisq for scat
#     - variance explained: R^2 (gaussian) or deviance explained (scat)
#     - 1000-iteration bootstrap CIs (coefboot.pffr), set.seed(42)
#     - globalbrainscore2 square-root transformed when BATCH_apply_gba_sqrt
#   Beta(d) plots with bootstrap CIs are produced at the end of every run.
#
#   OUTPUT (inside BATCH_output_dir)
#     RUN_INFO.txt, RUN_COMPLETE.txt
#     beta_plots/    beta_<term>.png for every term, all_terms.png
#     tables/        beta_curves.csv, beta_significance.csv, predictor_tests.csv,
#                    summary.csv, outliers_flagged.csv
#     exploratory/   raw_curves.png, mean_sd.png, functional_boxplot.png
#     diagnostics/   pffr_diagnostics.png, residuals.png, R2_by_density.png,
#                    pffr_coefficients.png
#     model/         FDA_results.RData
#
#   COMPUTATIONAL COST
#     Gaussian-family runs take several hours each on a 48-core HPC node,
#     including the 1000-iteration bootstrap; small-worldness under the scat
#     family can take up to a week per run. The saved outputs in results/ are
#     the recommended starting point for regenerating tables and figures.
#
################################################################################

if (!(exists("BATCH_MODE") && BATCH_MODE)) {
  stop("Run this through TWO_fda_group_difference_runs.R or THREE_fda_univariate_runs.R")
}

# ==============================================================================
# CONFIGURATION (from the runs driver)
# ==============================================================================
metric               <- BATCH_metric
all_predictors       <- BATCH_all_predictors
vars_to_skip_scaling <- BATCH_vars_to_skip_scaling
density_min          <- BATCH_density_min
density_max          <- BATCH_density_max
n_bootstrap          <- BATCH_n_bootstrap
p_threshold          <- BATCH_p_threshold
pffr_k_basis         <- BATCH_pffr_k_basis
pffr_family          <- BATCH_pffr_family
subjects_to_exclude  <- BATCH_subjects_to_exclude
apply_gba_sqrt       <- BATCH_apply_gba_sqrt
outlier_percentile   <- BATCH_outlier_percentile
output_dir           <- BATCH_output_dir
data_path            <- BATCH_data_path
run_meta             <- BATCH_run_meta          # list: pipeline, analysis, metric_label, model_label, index
predictor_labels     <- BATCH_predictor_labels  # named character vector for plot titles

n_cores <- suppressWarnings(as.integer(Sys.getenv("LSB_DJOB_NUMPROC")))
if (is.na(n_cores) || n_cores < 1) n_cores <- max(1, parallel::detectCores() - 1)

# ==============================================================================
# SETUP
# ==============================================================================
sub_dirs <- c("beta_plots", "tables", "exploratory", "diagnostics", "model")
for (d in c(output_dir, file.path(output_dir, sub_dirs))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
out <- function(sub, name) file.path(output_dir, sub, name)
unlink(file.path(output_dir, "RUN_COMPLETE.txt"))   # stale marker from an earlier attempt

suppressPackageStartupMessages({
  library(refund);    library(fda);     library(mgcv);      library(readxl)
  library(tidyverse); library(ggplot2); library(viridis);   library(gridExtra)
  library(fda.usc);   library(parallel)
})

set.seed(42)
use_F_test <- inherits(pffr_family, "family") && pffr_family$family == "gaussian"
family_name <- if (use_F_test) "gaussian" else "scat"
plab <- function(v) if (v %in% names(predictor_labels)) unname(predictor_labels[[v]]) else v
run_start <- Sys.time()

cat("\n==========================================\n")
cat("RUN CONFIGURATION\n")
cat("==========================================\n")
cat(sprintf("Index:        %s\n", run_meta$index))
cat(sprintf("Pipeline:     %s\n", run_meta$pipeline))
cat(sprintf("Analysis:     %s\n", run_meta$analysis))
cat(sprintf("Metric:       %s (%s)\n", run_meta$metric_label, metric))
cat(sprintf("Model:        %s\n", run_meta$model_label))
cat(sprintf("Predictors:   %s\n", paste(all_predictors, collapse = ", ")))
cat(sprintf("Not scaled:   %s\n", paste(intersect(all_predictors, vars_to_skip_scaling), collapse = ", ")))
cat(sprintf("Family:       %s\n", family_name))
cat(sprintf("Density:      %g-%g%%\n", density_min, density_max))
cat(sprintf("Excluded IDs: %s\n", if (length(subjects_to_exclude)) paste(subjects_to_exclude, collapse = ", ") else "none"))
cat(sprintf("Bootstrap:    %d iterations, %d cores\n", n_bootstrap, n_cores))
cat(sprintf("Data:         %s\n", data_path))
cat(sprintf("Output:       %s\n\n", output_dir))

# ==============================================================================
# SECTION 1: DATA
# ==============================================================================
cat("==========================================\nSECTION 1: DATA\n==========================================\n\n")

df <- read_excel(data_path)
n_original <- nrow(df)

if (length(subjects_to_exclude) > 0) {
  df <- df[!df$ID %in% subjects_to_exclude, ]
  cat(sprintf("Excluded %d subject(s) by ID. N = %d\n", n_original - nrow(df), nrow(df)))
}

missing_pred <- setdiff(all_predictors, colnames(df))
if (length(missing_pred)) stop("Predictor column(s) not in data: ", paste(missing_pred, collapse = ", "))

gba_sqrt_applied <- FALSE
if (apply_gba_sqrt && "globalbrainscore2" %in% all_predictors) {
  cat(sprintf("globalbrainscore2 before sqrt: range [%.2f, %.2f]\n",
              min(df$globalbrainscore2, na.rm = TRUE), max(df$globalbrainscore2, na.rm = TRUE)))
  df$globalbrainscore2 <- sqrt(df$globalbrainscore2)
  gba_sqrt_applied <- TRUE
  cat(sprintf("globalbrainscore2 after sqrt:  range [%.2f, %.2f]\n",
              min(df$globalbrainscore2, na.rm = TRUE), max(df$globalbrainscore2, na.rm = TRUE)))
}

# --- functional response ---
metric_cols   <- grep(paste0("^", metric, "_[0-9]+\\.[0-9]+$"), colnames(df), value = TRUE)
if (length(metric_cols) == 0) stop("No columns matching ", metric, "_<density>")
densities_all <- as.numeric(sub(paste0("^", metric, "_"), "", metric_cols))
keep          <- densities_all >= density_min & densities_all <= density_max
densities     <- densities_all[keep]
metric_cols   <- metric_cols[keep]
Y <- as.matrix(df[, metric_cols]); rownames(Y) <- df$ID; colnames(Y) <- densities
cat(sprintf("Response: %s, %d densities (%.2f-%.2f%%), range [%.4f, %.4f]\n",
            metric, length(densities), min(densities), max(densities),
            min(Y, na.rm = TRUE), max(Y, na.rm = TRUE)))
if (any(is.na(Y))) stop(sprintf("%d missing values in the response within %g-%g%%",
                                sum(is.na(Y)), density_min, density_max))

# --- predictors ---
predictor_df <- df %>% dplyr::select(ID, all_of(all_predictors)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(.)))
cat("\nMissingness in predictors:\n")
for (v in all_predictors) {
  nm <- sum(is.na(predictor_df[[v]]))
  cat(sprintf("  %-22s %d\n", v, nm))
}
complete_cases        <- complete.cases(predictor_df)
Y_complete            <- Y[complete_cases, , drop = FALSE]
predictor_df_complete <- predictor_df[complete_cases, ]
cat(sprintf("Complete cases: %d / %d\n", nrow(Y_complete), nrow(predictor_df)))

continuous_vars     <- setdiff(all_predictors, vars_to_skip_scaling)
predictor_df_scaled <- predictor_df_complete
for (v in continuous_vars) predictor_df_scaled[[v]] <- scale(predictor_df_scaled[[v]])[, 1]
cat(sprintf("\nz-scored:   %s\n", if (length(continuous_vars)) paste(continuous_vars, collapse = ", ") else "(none)"))
cat(sprintf("Not scaled: %s\n", paste(intersect(all_predictors, vars_to_skip_scaling), collapse = ", ")))

Y_clean                   <- Y_complete
predictor_df_clean        <- predictor_df_complete
predictor_df_scaled_clean <- predictor_df_scaled
N_final <- nrow(Y_clean)
cat(sprintf("\nFinal analysis sample: %d subjects\n", N_final))

# --- run info file ---
writeLines(c(
  sprintf("Index:          %s", run_meta$index),
  sprintf("Pipeline:       %s", run_meta$pipeline),
  sprintf("Analysis:       %s", run_meta$analysis),
  sprintf("Metric:         %s (%s)", run_meta$metric_label, metric),
  sprintf("Model:          %s", run_meta$model_label),
  sprintf("Predictors:     %s", paste(all_predictors, collapse = ", ")),
  sprintf("z-scored:       %s", paste(continuous_vars, collapse = ", ")),
  sprintf("Not scaled:     %s", paste(intersect(all_predictors, vars_to_skip_scaling), collapse = ", ")),
  sprintf("GBA sqrt:       %s", gba_sqrt_applied),
  sprintf("Family:         %s", family_name),
  sprintf("Density range:  %g-%g%% (%d points)", density_min, density_max, length(densities)),
  sprintf("Excluded IDs:   %s", if (length(subjects_to_exclude)) paste(subjects_to_exclude, collapse = ", ") else "none"),
  sprintf("N:              %d (of %d in file)", N_final, n_original),
  sprintf("Bootstrap:      %d iterations", n_bootstrap),
  sprintf("Data file:      %s", data_path),
  sprintf("Started:        %s", format(run_start)),
  sprintf("Host:           %s", Sys.info()[["nodename"]])
), file.path(output_dir, "RUN_INFO.txt"))

# ==============================================================================
# SECTION 2: OUTLIER FLAGGING (functional depth; flagged only, nobody removed)
# ==============================================================================
cat("\n==========================================\nSECTION 2: OUTLIER FLAGGING (not removed)\n==========================================\n\n")
outlier_info <- data.frame()
depth_values <- NULL
tryCatch({
  depth_values <- depth.mode(fdata(Y_clean, argvals = densities))$dep
  thr <- quantile(depth_values, outlier_percentile)
  idx <- which(depth_values < thr)
  outlier_info <- data.frame(ID = predictor_df_clean$ID[idx], Depth = depth_values[idx],
                             Max_Value = apply(Y_clean[idx, , drop = FALSE], 1, max),
                             Min_Value = apply(Y_clean[idx, , drop = FALSE], 1, min),
                             Mean_Value = rowMeans(Y_clean[idx, , drop = FALSE]))
  cat(sprintf("Lowest %.0f%% depth (depth < %.4f): %d subject(s): %s\n", outlier_percentile * 100, thr,
              length(idx), paste(outlier_info$ID, collapse = ", ")))
  write.csv(outlier_info, out("tables", "outliers_flagged.csv"), row.names = FALSE)
}, error = function(e) cat("Depth computation failed:", conditionMessage(e), "\n"))

# ==============================================================================
# SECTION 3: EXPLORATORY PLOTS
# ==============================================================================
cat("\n==========================================\nSECTION 3: EXPLORATORY PLOTS\n==========================================\n\n")
title_tag <- sprintf("%s | %s | %s", run_meta$pipeline, run_meta$metric_label, run_meta$model_label)

try({
  plot_df <- as.data.frame(Y_clean) %>% mutate(ID = predictor_df_clean$ID) %>%
    pivot_longer(-ID, names_to = "density", values_to = "value") %>% mutate(density = as.numeric(density))
  p1 <- ggplot(plot_df, aes(density, value, group = ID)) + geom_line(alpha = 0.3, linewidth = 0.3) +
    labs(x = "Network Density (%)", y = run_meta$metric_label,
         title = sprintf("Raw %s curves (n = %d)", run_meta$metric_label, N_final), subtitle = title_tag) +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  ggsave(out("exploratory", "raw_curves.png"), p1, width = 10, height = 6, dpi = 300)
  
  mean_df <- data.frame(density = densities, mean = colMeans(Y_clean), sd = apply(Y_clean, 2, sd))
  p2 <- ggplot(mean_df, aes(density)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "steelblue", alpha = 0.3) +
    geom_line(aes(y = mean), color = "darkblue", linewidth = 1) +
    labs(x = "Network Density (%)", y = run_meta$metric_label,
         title = sprintf("Mean ± SD %s", run_meta$metric_label), subtitle = title_tag) +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  ggsave(out("exploratory", "mean_sd.png"), p2, width = 10, height = 6, dpi = 300)
  cat("  raw_curves.png, mean_sd.png\n")
})

try({
  png(out("exploratory", "functional_boxplot.png"), width = 10, height = 6, units = "in", res = 300)
  par(mar = c(5, 4, 4, 2))
  fbplot(t(Y_clean), x = densities, xlab = "Network Density (%)", ylab = run_meta$metric_label,
         main = sprintf("Functional boxplot: %s", title_tag), color = viridis(3)[2], xlim = range(densities))
  dev.off()
  cat("  functional_boxplot.png\n")
})
if (dev.cur() > 1) dev.off()

# ==============================================================================
# SECTION 4: PFFR MODEL
# ==============================================================================
cat("\n==========================================\nSECTION 4: PFFR MODEL\n==========================================\n\n")
pffr_data   <- predictor_df_scaled_clean
pffr_data$Y <- Y_clean
bsy <- list(bs = "ps", k = pffr_k_basis, m = c(2, 1))

formula_parts <- paste(all_predictors, collapse = " + ")
pffr_formula  <- as.formula(paste("Y ~", formula_parts))
cat(sprintf("Formula: Y ~ %s | family: %s\n\n", formula_parts, family_name))

pffr_fit <- pffr(pffr_formula, yind = densities, data = pffr_data, bs.yindex = bsy, family = pffr_family)
print(summary(pffr_fit))

# ==============================================================================
# SECTION 5: HYPOTHESIS TESTS
# ==============================================================================
cat("\n==========================================\nSECTION 5: HYPOTHESIS TESTS\n==========================================\n\n")

get_var_explained <- function(model) { s <- summary(model); if (use_F_test) s$r.sq else s$dev.expl }
ve_label <- if (use_F_test) "R^2" else "Deviance explained"

pffr_intercept_only <- pffr(Y ~ 1, yind = densities, data = pffr_data, bs.yindex = bsy, family = pffr_family)
ve_full       <- get_var_explained(pffr_fit)
ve_intercept  <- get_var_explained(pffr_intercept_only)
functional_ve <- (ve_full - ve_intercept) / (1 - ve_intercept)
cat(sprintf("%s full model: %.4f | functional %s: %.4f\n\n", ve_label, ve_full, ve_label, functional_ve))

# overall LR test: varying vs constant coefficients
constant_formula <- as.formula(paste("Y ~ 1 +", paste0("c(", all_predictors, ")", collapse = " + ")))
pffr_constant <- pffr(constant_formula, yind = densities, data = pffr_data, bs.yindex = bsy, family = pffr_family)
mtest_full <- pffr_fit; mtest_constant <- pffr_constant
class(mtest_full) <- class(mtest_full)[-1]; class(mtest_constant) <- class(mtest_constant)[-1]
lr_test_overall <- anova(mtest_constant, mtest_full, test = if (use_F_test) "F" else "Chisq")
overall_stat <- if (use_F_test) lr_test_overall$F[2] else lr_test_overall$Deviance[2]
overall_p    <- if (use_F_test) lr_test_overall$`Pr(>F)`[2] else lr_test_overall$`Pr(>Chi)`[2]
overall_model_significant <- overall_p < p_threshold
cat(sprintf("Overall LR test: %s = %.4f, p = %.3e, significant at p < %.3f: %s\n\n",
            if (use_F_test) "F" else "Chisq", overall_stat, overall_p, p_threshold,
            ifelse(overall_model_significant, "YES", "NO")))

# per-predictor LR tests
predictor_tests <- data.frame(Predictor = character(), Label = character(), Test_Statistic = numeric(),
                              Test_Type = character(), p_value = numeric(), functional_ve = numeric(),
                              semi_partial_r = numeric(), stringsAsFactors = FALSE)
if (length(all_predictors) == 1) {
  predictor_tests <- data.frame(Predictor = all_predictors, Label = plab(all_predictors),
                                Test_Statistic = overall_stat, Test_Type = if (use_F_test) "F" else "Chisq",
                                p_value = overall_p, functional_ve = functional_ve,
                                semi_partial_r = sqrt(max(functional_ve, 0)), stringsAsFactors = FALSE)
  cat("Single predictor: per-predictor test = overall test\n")
} else {
  for (pred in all_predictors) {
    red <- tryCatch(pffr(as.formula(paste("Y ~", paste(setdiff(all_predictors, pred), collapse = " + "))),
                         yind = densities, data = pffr_data, bs.yindex = bsy, family = pffr_family),
                    error = function(e) NULL)
    if (is.null(red)) { cat(sprintf("  %s: reduced model failed\n", pred)); next }
    mr <- red; class(mr) <- class(mr)[-1]
    lr <- tryCatch(anova(mr, mtest_full, test = if (use_F_test) "F" else "Chisq"), error = function(e) NULL)
    if (is.null(lr)) { cat(sprintf("  %s: LR test failed\n", pred)); next }
    single <- tryCatch(pffr(as.formula(paste("Y ~ 1 +", pred)), yind = densities, data = pffr_data,
                            bs.yindex = bsy, family = pffr_family), error = function(e) NULL)
    fve <- if (!is.null(single)) (get_var_explained(single) - ve_intercept) / (1 - ve_intercept) else NA
    sr  <- sqrt(abs(ve_full - get_var_explained(red)))
    st  <- if (use_F_test) lr$F[2] else lr$Deviance[2]
    pv  <- if (use_F_test) lr$`Pr(>F)`[2] else lr$`Pr(>Chi)`[2]
    predictor_tests <- rbind(predictor_tests, data.frame(Predictor = pred, Label = plab(pred), Test_Statistic = st,
                                                         Test_Type = if (use_F_test) "F" else "Chisq", p_value = pv,
                                                         functional_ve = fve, semi_partial_r = sr,
                                                         stringsAsFactors = FALSE))
    cat(sprintf("  %-24s %s = %9.3f, p = %.3e, functional %s = %.4f, sr = %.4f%s\n", pred,
                if (use_F_test) "F" else "Chisq", st, pv, ve_label, fve, sr,
                ifelse(pv < p_threshold, " ***", ifelse(pv < 0.01, " **", ifelse(pv < 0.05, " *", "")))))
  }
}
write.csv(predictor_tests, out("tables", "predictor_tests.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 6: BOOTSTRAP CIs
# ==============================================================================
cat("\n==========================================\nSECTION 6: BOOTSTRAP CIs\n==========================================\n\n")
cat(sprintf("Running %d bootstrap iterations on %d cores (started %s)\n", n_bootstrap, n_cores, format(Sys.time())))
bootstrap_coefs <- tryCatch(coefboot.pffr(pffr_fit, B = n_bootstrap, ncpus = n_cores, parallel = "multicore"),
                            error = function(e) { cat("coefboot.pffr error:", conditionMessage(e), "\n"); NULL })
cat(sprintf("Bootstrap %s (%s)\n", if (is.null(bootstrap_coefs)) "FAILED" else "complete", format(Sys.time())))

# ==============================================================================
# SECTION 7: BETA(d) PLOTS AND TABLES
# ==============================================================================
cat("\n==========================================\nSECTION 7: BETA(d) PLOTS\n==========================================\n\n")

getPlotObject <- function(model) { ff <- tempfile(); svg(filename = ff); po <- plot(model); dev.off(); unlink(ff); po }

fmt_d <- function(d) if (abs(d - round(d)) < 1e-6) sprintf("%d", as.integer(round(d))) else sprintf("%.1f", d)
sig_ranges <- function(d, lo, hi) {
  ex <- (lo > 0) | (hi < 0)
  if (!any(ex, na.rm = TRUE)) return("ns")
  r <- rle(ex); ends <- cumsum(r$lengths); starts <- c(1L, head(ends, -1) + 1L)
  segs <- mapply(function(s, e, v) if (isTRUE(v)) sprintf("%s-%s%%", fmt_d(d[s]), fmt_d(d[e])) else NA,
                 starts, ends, r$values)
  paste(na.omit(segs), collapse = ", ")
}

term_names <- c("Intercept", all_predictors)
beta_curves <- data.frame(); beta_sig <- data.frame()

if (!is.null(bootstrap_coefs)) tryCatch({
  po  <- getPlotObject(pffr_fit)
  sm  <- bootstrap_coefs$smterms
  if (length(po) != length(sm) || length(sm) != length(term_names)) {
    cat(sprintf("WARNING: %d plot terms, %d bootstrap terms, %d expected terms; bootstrap term names: %s\n",
                length(po), length(sm), length(term_names), paste(names(sm), collapse = ", ")))
  }
  plots <- list()
  for (i in seq_len(min(length(po), length(sm), length(term_names)))) {
    tn <- term_names[i]
    # point estimate from the fitted model (as in the paper's figures); if the plot
    # grid and bootstrap grid differ in length, fall back to the bootstrap object's
    # own estimate and grid so the run still produces output
    if (length(po[[i]]$x) == length(sm[[i]][["2.5%"]])) {
      bd <- data.frame(density = as.numeric(po[[i]]$x), beta = as.numeric(po[[i]]$fit),
                       ci_lower = sm[[i]][["2.5%"]], ci_upper = sm[[i]][["97.5%"]])
    } else {
      cat(sprintf("  NOTE: %s plot grid (%d) != bootstrap grid (%d); using bootstrap estimate\n",
                  tn, length(po[[i]]$x), length(sm[[i]][["2.5%"]])))
      bd <- data.frame(density = sm[[i]][[2]], beta = sm[[i]][[1]],
                       ci_lower = sm[[i]][["2.5%"]], ci_upper = sm[[i]][["97.5%"]])
    }
    bd$ci_excludes_zero <- (bd$ci_lower > 0) | (bd$ci_upper < 0)
    rng <- sig_ranges(bd$density, bd$ci_lower, bd$ci_upper)
    beta_curves <- rbind(beta_curves, data.frame(term = tn, label = if (tn == "Intercept") "Intercept" else plab(tn),
                                                 bootstrap_term = names(sm)[i], bd))
    sig <- bd$ci_excludes_zero
    beta_sig <- rbind(beta_sig, data.frame(
      term = tn, label = if (tn == "Intercept") "Intercept" else plab(tn), bootstrap_term = names(sm)[i],
      significant = any(sig), ci_excludes_zero_ranges = rng,
      pct_density_significant = round(100 * mean(sig), 1),
      direction_in_sig_region = if (any(sig)) ifelse(mean(bd$beta[sig]) > 0, "positive", "negative") else NA,
      mean_beta = mean(bd$beta), stringsAsFactors = FALSE))
    
    lab <- if (tn == "Intercept") "Intercept" else plab(tn)
    gg <- ggplot(bd, aes(density, beta)) +
      geom_hline(yintercept = 0, lty = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") +
      geom_line(linewidth = 1, color = "darkblue") +
      ylab(expression(beta(d))) + xlab("Density (%)") +
      ggtitle(sprintf("%s effect on %s", lab, run_meta$metric_label),
              subtitle = sprintf("%s | %s | N = %d | CI excludes 0: %s", run_meta$pipeline,
                                 run_meta$model_label, N_final, rng)) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 8))
    fn <- sprintf("beta_%s.png", gsub("[^A-Za-z0-9]+", "_", tn))
    ggsave(out("beta_plots", fn), gg, width = 8, height = 5, dpi = 300)
    plots[[length(plots) + 1]] <- gg
    cat(sprintf("  %-26s CI excludes 0: %s\n", fn, rng))
  }
  nc <- ceiling(sqrt(length(plots))); nr <- ceiling(length(plots) / nc)
  ggsave(out("beta_plots", "all_terms.png"),
         arrangeGrob(grobs = plots, ncol = nc, top = title_tag),
         width = 7 * nc, height = 4.5 * nr, dpi = 200, limitsize = FALSE)
  cat("  all_terms.png\n")
  write.csv(beta_curves, out("tables", "beta_curves.csv"), row.names = FALSE)
  write.csv(beta_sig, out("tables", "beta_significance.csv"), row.names = FALSE)
}, error = function(e) cat("Beta plotting failed:", conditionMessage(e), "\n"))
if (is.null(bootstrap_coefs)) cat("No bootstrap: beta plots and tables skipped.\n")

# ==============================================================================
# SECTION 8: DIAGNOSTICS
# ==============================================================================
cat("\n==========================================\nSECTION 8: DIAGNOSTICS\n==========================================\n\n")

try({
  png(out("diagnostics", "pffr_coefficients.png"), width = 14, height = 12, units = "in", res = 300)
  n_plots <- length(all_predictors) + 1; nc <- ceiling(sqrt(n_plots)); nr <- ceiling(n_plots / nc)
  par(mfrow = c(nr, nc))
  plot(pffr_fit, pages = 1, shade = TRUE, shade.col = "grey80", all.terms = TRUE, rug = FALSE)
  dev.off(); cat("  pffr_coefficients.png\n")
})
if (dev.cur() > 1) dev.off()

try({
  png(out("diagnostics", "pffr_diagnostics.png"), width = 12, height = 10, units = "in", res = 300)
  pffr.check(pffr_fit)
  dev.off(); cat("  pffr_diagnostics.png\n")
})
if (dev.cur() > 1) dev.off()

R2_function <- NULL
try({
  residuals_pffr <- Y_clean - fitted(pffr_fit)
  resid_df <- as.data.frame(residuals_pffr) %>% mutate(ID = predictor_df_clean$ID) %>%
    pivot_longer(-ID, names_to = "density", values_to = "residual") %>% mutate(density = as.numeric(density))
  pr <- ggplot(resid_df, aes(density, residual, group = ID)) + geom_line(alpha = 0.2) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(x = "Network Density (%)", y = "Residual (observed - fitted)", title = "Residual curves", subtitle = title_tag) +
    theme_minimal()
  ggsave(out("diagnostics", "residuals.png"), pr, width = 10, height = 6, dpi = 300)
  R2_function <- 1 - colSums(residuals_pffr^2) / colSums(scale(Y_clean, scale = FALSE)^2)
  pR <- ggplot(data.frame(density = densities, R2 = R2_function), aes(density, R2)) +
    geom_line(color = "darkgreen", linewidth = 1) + geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Network Density (%)", y = expression(R^2), title = "R² by density", subtitle = title_tag) +
    theme_minimal()
  ggsave(out("diagnostics", "R2_by_density.png"), pR, width = 10, height = 6, dpi = 300)
  cat(sprintf("  residuals.png, R2_by_density.png (mean R^2 across density %.4f)\n", mean(R2_function)))
})

# ==============================================================================
# SECTION 9: SAVE
# ==============================================================================
cat("\n==========================================\nSECTION 9: SAVE\n==========================================\n\n")

summary_results <- data.frame(
  Index = run_meta$index, Pipeline = run_meta$pipeline, Analysis = run_meta$analysis,
  Metric = metric, Metric_Label = run_meta$metric_label, Model = run_meta$model_label,
  Predictors = paste(all_predictors, collapse = ", "), Family = family_name,
  Density_Range = paste0(density_min, "-", density_max, "%"), N_Original = n_original,
  N_Excluded_IDs = length(subjects_to_exclude), N_Final = N_final,
  N_Depth_Flagged = nrow(outlier_info), GBA_Sqrt_Applied = gba_sqrt_applied,
  Variance_Explained = ve_full, Functional_VE = functional_ve, VE_Label = ve_label,
  Overall_Test_Statistic = overall_stat, Overall_Test_Type = if (use_F_test) "F" else "Chisq",
  Overall_p = overall_p, Overall_Significant = overall_model_significant, P_Threshold = p_threshold,
  Bootstrap_Completed = !is.null(bootstrap_coefs), stringsAsFactors = FALSE)
write.csv(summary_results, out("tables", "summary.csv"), row.names = FALSE)

save(metric, densities, density_min, density_max, pffr_family, family_name,
     Y_clean, Y_complete, predictor_df_clean, predictor_df_scaled_clean,
     pffr_fit, pffr_intercept_only, pffr_constant, depth_values, outlier_info,
     predictor_tests, summary_results, lr_test_overall, ve_full, ve_intercept, functional_ve, ve_label,
     bootstrap_coefs, beta_curves, beta_sig, R2_function,
     subjects_to_exclude, all_predictors, vars_to_skip_scaling, gba_sqrt_applied, run_meta,
     file = out("model", "FDA_results.RData"))
cat("  model/FDA_results.RData\n")

writeLines(c(sprintf("Finished: %s", format(Sys.time())),
             sprintf("Elapsed:  %.1f hours", as.numeric(difftime(Sys.time(), run_start, units = "hours"))),
             sprintf("Bootstrap completed: %s", !is.null(bootstrap_coefs))),
           file.path(output_dir, "RUN_COMPLETE.txt"))

cat("\n==========================================\nRUN SUMMARY\n==========================================\n")
cat(sprintf("%s | %s | %s | %s | N = %d | family %s\n", run_meta$pipeline, run_meta$analysis,
            run_meta$metric_label, run_meta$model_label, N_final, family_name))
cat(sprintf("Overall: %s = %.4f, p = %.3e (%s)\n", if (use_F_test) "F" else "Chisq", overall_stat, overall_p,
            ifelse(overall_model_significant, "significant", "not significant")))
if (nrow(beta_sig)) print(beta_sig[, c("term", "significant", "ci_excludes_zero_ranges", "direction_in_sig_region")], row.names = FALSE)
cat("\nSession info:\n"); print(sessionInfo())