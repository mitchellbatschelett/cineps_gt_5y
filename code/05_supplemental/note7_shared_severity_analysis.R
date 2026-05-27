################################################################################
#
#   SUPPLEMENTARY NOTE 7: SHARED-SEVERITY ASSESSMENT ANALYSIS
#
#   Standalone analysis generator for Supplementary Note 7 (Figs 17-20, Tbl 7).
#   Generates the Supplementary Note 7 outputs from the analytic cohort:
#     - repo-relative paths (auto-detected repo_root)
#     - data source = data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx
#     - all outputs written to results/post_lasso_shared_severity/
#     - ADDS a per-subject PC1 + FPC1 + FPC2 score export (PC1_FPC_scores.csv)
#       so the supplement builder can render Supp Fig 19 by reading scores
#       instead of re-running the FPCA pipeline.
#
#   OUTPUTS (results/post_lasso_shared_severity/)
#     exposure_correlation_matrix.csv          (Supp Fig 17)
#     exposure_clusters.csv                     (context; Fig 17 prose)
#     exposure_VIFs.csv, forced_covariate_r2.csv (context)
#     exposure_pca_variance.csv                 (Supp Fig 18A)
#     exposure_pca_loadings.csv                 (Supp Fig 18B)
#     PC1_FPC_scores.csv                        (Supp Fig 19; NEW)
#     perm_global_exposure_set_by_metric.csv    (Supp Fig 20A, Tbl 7)
#     perm_incremental_exposure_by_metric.csv   (Supp Fig 20B, Tbl 7)
#     perm_latent_severity_by_metric.csv        (Supp Fig 20C, Tbl 7)
#
#   USAGE
#     Rscript code/05_supplemental/note7_shared_severity_analysis.R
#     (or source() from an R session with repo_root pre-set)
#
#   REQUIRES  R 4.4.0; refund, mgcv, readxl, tidyverse, car, fda.usc
#
################################################################################

# ==============================================================================
# REPO ROOT (auto-detect)
# ==============================================================================

.is_repo_root <- function(p) {
  dir.exists(file.path(p, "data", "analysis_ready")) &&
    dir.exists(file.path(p, "code"))
}
.find_repo_root_from_path <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  for (i in 1:8) {
    if (.is_repo_root(p)) return(p)
    parent <- dirname(p); if (parent == p) break; p <- parent
  }
  NULL
}
if (exists("repo_root", inherits = TRUE) &&
    is.character(repo_root) && length(repo_root) == 1 &&
    .is_repo_root(normalizePath(repo_root, winslash = "/", mustWork = FALSE))) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
} else {
  repo_root <- NULL
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    this_script <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
    repo_root <- .find_repo_root_from_path(dirname(this_script))
  }
  if (is.null(repo_root)) repo_root <- .find_repo_root_from_path(getwd())
  if (is.null(repo_root))
    stop("Could not locate repo root. Set it manually before running:\n",
         "  repo_root <- \"/full/path/to/repository\"")
}
cat(sprintf("Repo root: %s\n\n", repo_root))

# ==============================================================================
# CONFIGURATION
# ==============================================================================

data_path  <- file.path(repo_root, "data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx")
output_dir <- file.path(repo_root, "results/post_lasso_shared_severity")

forced_covariates <- c("eTIV", "sex", "sriskscore")
exposures <- c("bpd2", "bw_z", "ga", "globalbrainscore2", "anyrop", "sepsis2")

categorical_exposures <- c("bpd2", "anyrop", "sepsis2")
categorical_forced    <- c("sex")

gbs2_transform <- "sqrt"   # "sqrt", "log1p", or "none"

metrics <- c("rand_norm_wei_ACC", "rand_norm_wei_GE", "str", "rand_norm_wei_SW")
metric_labels <- c(
  "rand_norm_wei_ACC" = "Clustering (ACC)",
  "rand_norm_wei_GE"  = "Global Efficiency (GE)",
  "str"               = "Strength",
  "rand_norm_wei_SW"  = "Small-Worldness (SW)"
)

metric_subjects_to_exclude <- list(
  "rand_norm_wei_ACC" = c(),
  "rand_norm_wei_GE"  = c(659),
  "str"               = c(128, 713),
  "rand_norm_wei_SW"  = c(309, 321, 8155)
)

exposure_labels <- c(
  "bpd2" = "BPD",
  "bw_z" = "Birth Weight Z",
  "ga" = "Gestational Age",
  "globalbrainscore2" = "GBS2",
  "anyrop" = "ROP",
  "sepsis2" = "Sepsis"
)

density_min <- 11
density_max <- 100

fpca_pve_threshold <- 0.995
fpca_min_npc <- 3
fpca_max_npc <- 10

pffr_k_basis <- 20
pffr_family  <- gaussian()

remove_outliers    <- FALSE
outlier_percentile <- 0.01

B_perm <- 1000
set.seed(42)

cluster_method <- "hclust"
cluster_dist   <- "1-abs(cor)"
cluster_cut    <- "k"
cluster_k      <- 2
cluster_height <- 0.55
cluster_min_abs_cor_for_report <- 0.40

use_selection_diagnostics <- FALSE   # selection diagnostics not part of Note 7 outputs

# ==============================================================================
# LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(refund)   # fpca.sc + pffr
  library(mgcv)     # backend for pffr
  library(readxl)
  library(tidyverse)
  library(car)      # vif
  library(fda.usc)  # functional outlier depth (only if remove_outliers)
})

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# HELPERS  (preserved from source)
# ==============================================================================

get_density_cols <- function(df, metric_prefix, d_min, d_max) {
  pat <- paste0("^", metric_prefix, "_[0-9]+(\\.[0-9]+)?$")
  cands <- grep(pat, names(df), value = TRUE)
  if (length(cands) == 0) return(character(0))
  d_vals <- as.numeric(sub(paste0("^", metric_prefix, "_"), "", cands))
  keep <- !is.na(d_vals) & d_vals >= d_min & d_vals <= d_max
  cands[keep][order(d_vals[keep])]
}

# Freedman-Lane permutation test for incremental contribution of a predictor
# set; compares full vs reduced design in multivariate response S (n x K).
perm_test_increment <- function(S, X_full, X_red, B = 1000) {
  XtX_full <- crossprod(X_full)
  XtX_red  <- crossprod(X_red)
  H_full <- X_full %*% solve(XtX_full) %*% t(X_full)
  H_red  <- X_red  %*% solve(XtX_red ) %*% t(X_red )
  fit_red <- lm.fit(x = X_red, y = S)
  E <- fit_red$residuals
  S_hat_red <- X_red %*% coef(fit_red)
  T_obs <- sum(diag(t(S) %*% (H_full - H_red) %*% S))
  T_perm <- replicate(B, {
    idx <- sample.int(nrow(S))
    S_star <- S_hat_red + E[idx, , drop = FALSE]
    sum(diag(t(S_star) %*% (H_full - H_red) %*% S_star))
  })
  p <- (1 + sum(T_perm >= T_obs)) / (B + 1)
  list(T_obs = T_obs, p_perm = p)
}

# ==============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ==============================================================================

cat("--- Loading Data ---\n")
df <- as.data.frame(read_excel(data_path))

if (gbs2_transform == "sqrt") {
  df$globalbrainscore2 <- sqrt(df$globalbrainscore2)
  cat("Applied sqrt transform to globalbrainscore2\n")
} else if (gbs2_transform == "log1p") {
  df$globalbrainscore2 <- log1p(df$globalbrainscore2)
  cat("Applied log1p transform to globalbrainscore2\n")
} else {
  cat("No transform applied to globalbrainscore2\n")
}

all_vars <- c(forced_covariates, exposures)
complete_mask <- complete.cases(df[, all_vars])
df_complete <- df[complete_mask, ]
cat(sprintf("Complete cases (forced + exposures): %d / %d\n\n",
            nrow(df_complete), nrow(df)))

df_complete <- df_complete %>%
  mutate(across(all_of(categorical_exposures), ~ as.integer(.x))) %>%
  mutate(across(all_of(categorical_forced), ~ factor(.x)))

df_std <- df_complete
for (v in c(forced_covariates, exposures)) {
  if (!(v %in% c(categorical_exposures, categorical_forced))) {
    df_std[[v]] <- as.numeric(scale(df_std[[v]])[, 1])
  }
}

# ==============================================================================
# SECTION 2: EXPOSURE INTERCORRELATION + CLUSTERING (ONCE)
# ==============================================================================

cat("========================================\n")
cat("SECTION 2: EXPOSURE INTERCORRELATION + CLUSTERING\n")
cat("========================================\n\n")

exp_mat <- as.matrix(df_std[, exposures, drop = FALSE])
colnames(exp_mat) <- exposure_labels[exposures]
corr_mat <- cor(exp_mat, use = "pairwise.complete.obs")
print(round(corr_mat, 3))

cat("\nPairwise |r| >", cluster_min_abs_cor_for_report, ":\n")
for (i in 1:(ncol(corr_mat) - 1)) {
  for (j in (i + 1):ncol(corr_mat)) {
    if (abs(corr_mat[i, j]) > cluster_min_abs_cor_for_report) {
      cat(sprintf("  %-18s x %-18s  r = %.3f\n",
                  colnames(corr_mat)[i], colnames(corr_mat)[j], corr_mat[i, j]))
    }
  }
}

write.csv(as.data.frame(round(corr_mat, 4)),
          file.path(output_dir, "exposure_correlation_matrix.csv"))
cat("\nSaved: exposure_correlation_matrix.csv\n")

dist_mat <- as.dist(1 - abs(corr_mat))
hc <- hclust(dist_mat, method = "average")
if (cluster_cut == "k") {
  cl <- cutree(hc, k = cluster_k)
  cat(sprintf("\nExposure clusters via hclust cutree(k=%d):\n", cluster_k))
} else {
  cl <- cutree(hc, h = cluster_height)
  cat(sprintf("\nExposure clusters via hclust cutree(h=%.2f):\n", cluster_height))
}
labels_order <- exposure_labels[exposures]
cl_by_label <- cl[labels_order]
cluster_df <- tibble(exposure = exposures, label = labels_order,
                     cluster = as.integer(cl_by_label))
print(cluster_df)
write.csv(cluster_df, file.path(output_dir, "exposure_clusters.csv"), row.names = FALSE)
cat("Saved: exposure_clusters.csv\n")

# ==============================================================================
# SECTION 3: VIFs + Forced-covariate R^2 (context only)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 3: VIFs + Forced-covariate R^2 \n")
cat("========================================\n\n")

vif_formula <- as.formula(paste("rnorm(nrow(df_std)) ~",
                                paste(exposures, collapse = " + ")))
vif_model <- lm(vif_formula, data = df_std)
vif_vals <- car::vif(vif_model)
vif_df <- tibble(Exposure = exposure_labels[names(vif_vals)],
                 VIF = round(as.numeric(vif_vals), 2)) %>% arrange(desc(VIF))
print(vif_df)
write.csv(vif_df, file.path(output_dir, "exposure_VIFs.csv"), row.names = FALSE)
cat("Saved: exposure_VIFs.csv\n")

forced_r2 <- map_dfr(exposures, function(exp) {
  fmla <- as.formula(paste(exp, "~", paste(forced_covariates, collapse = " + ")))
  fit <- lm(fmla, data = df_std)
  tibble(Exposure = exposure_labels[exp], R2 = summary(fit)$r.squared)
}) %>% mutate(R2 = round(R2, 4)) %>% arrange(desc(R2))
print(forced_r2)
write.csv(forced_r2, file.path(output_dir, "forced_covariate_r2.csv"), row.names = FALSE)
cat("Saved: forced_covariate_r2.csv\n")

# ==============================================================================
# SECTION 4: PCA ON EXPOSURES (LATENT SEVERITY AXIS)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 4: PCA ON EXPOSURES (LATENT SEVERITY)\n")
cat("========================================\n\n")

exp_for_pca <- scale(df_complete[, exposures])
colnames(exp_for_pca) <- exposure_labels[exposures]
pca_exp <- prcomp(exp_for_pca, center = FALSE, scale. = FALSE)

pve <- pca_exp$sdev^2 / sum(pca_exp$sdev^2)
cum_pve <- cumsum(pve)

pca_results <- tibble(
  PC = paste0("PC", seq_along(pve)),
  Variance_Explained = round(pve * 100, 2),
  Cumulative = round(cum_pve * 100, 2)
)
write.csv(pca_results, file.path(output_dir, "exposure_pca_variance.csv"), row.names = FALSE)

loadings_df <- as.data.frame(round(pca_exp$rotation, 4)) %>%
  rownames_to_column("Exposure")
write.csv(loadings_df, file.path(output_dir, "exposure_pca_loadings.csv"), row.names = FALSE)
cat("Saved: exposure_pca_variance.csv, exposure_pca_loadings.csv\n")
cat(sprintf("  PC1 PVE = %.1f%% (n=%d)\n", 100 * pve[1], nrow(df_complete)))

# Subject-level PC1 / PC2 (full complete-case sample), for Supp Fig 19.
pc_scores_df <- data.frame(
  ID  = df_complete$ID,
  PC1 = as.numeric(exp_for_pca %*% pca_exp$rotation[, 1]),
  PC2 = as.numeric(exp_for_pca %*% pca_exp$rotation[, 2]),
  stringsAsFactors = FALSE
)

# ==============================================================================
# SECTION 5: PER-METRIC PFFR RESIDUALS -> FPCA -> FREEDMAN-LANE PERMUTATION
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 5: UNIFIED PER-METRIC ANALYSIS (FPCA + PERMUTATION)\n")
cat("========================================\n\n")

perm_results_all <- list()
global_results   <- list()
latent_results   <- list()
score_rows       <- list()   # NEW: per-subject PC1 + FPC scores for Fig 19

for (m in metrics) {
  
  cat(sprintf("\n--- %s (%s) ---\n", m, metric_labels[m]))
  
  d_cols <- get_density_cols(df_complete, m, density_min, density_max)
  if (length(d_cols) == 0) { cat("  No density columns found. Skipping.\n"); next }
  
  metric_mask <- complete.cases(df_complete[, d_cols])
  stopifnot(nrow(df_std) == nrow(df_complete))
  
  df_m <- df_std[metric_mask, , drop = FALSE]
  Y <- as.matrix(df_complete[metric_mask, d_cols, drop = FALSE])
  
  dens <- as.numeric(sub(paste0("^", m, "_"), "", d_cols))
  col_order <- order(dens); dens <- dens[col_order]
  Y <- Y[, col_order, drop = FALSE]; colnames(Y) <- dens
  
  excl_ids <- metric_subjects_to_exclude[[m]]
  if (!is.null(excl_ids) && length(excl_ids) > 0) {
    keep <- !df_m$ID %in% excl_ids
    n_excl <- sum(!keep)
    df_m <- df_m[keep, , drop = FALSE]; Y <- Y[keep, , drop = FALSE]
    cat(sprintf("  Excluded %d subjects by ID (%s). N = %d\n",
                n_excl, paste(excl_ids, collapse = ", "), nrow(df_m)))
  }
  
  if (remove_outliers) {
    fdata_obj <- fda.usc::fdata(Y, argvals = dens)
    depth_vals <- fda.usc::depth.mode(fdata_obj)$dep
    outlier_thresh <- quantile(depth_vals, outlier_percentile)
    outlier_idx <- which(depth_vals < outlier_thresh)
    if (length(outlier_idx) > 0) {
      Y <- Y[-outlier_idx, , drop = FALSE]
      df_m <- df_m[-outlier_idx, , drop = FALSE]
      cat(sprintf("  Removed %d functional outliers (depth < %.4f)\n",
                  length(outlier_idx), outlier_thresh))
    }
  }
  
  cat(sprintf("  N = %d, T = %d density points\n", nrow(df_m), ncol(Y)))
  
  # Forced-adjust curves using pffr
  pffr_data <- df_m; pffr_data$Y <- Y
  forced_formula <- as.formula(paste("Y ~", paste(forced_covariates, collapse = " + ")))
  pffr_fit <- pffr(
    forced_formula, yind = dens, data = pffr_data,
    bs.yindex = list(bs = "ps", k = pffr_k_basis, m = c(2, 1)),
    family = pffr_family
  )
  Y_resid <- residuals(pffr_fit); colnames(Y_resid) <- dens
  
  # FPCA on residual curves; truncate PCs (min/max + PVE threshold) like main.
  fp <- refund::fpca.sc(Y = Y_resid, argvals = dens,
                        pve = fpca_pve_threshold, npc = fpca_max_npc, var = TRUE)
  pve_ind <- fp$evalues / sum(fp$evalues)
  pve_cum <- cumsum(pve_ind)
  k_pve <- which(pve_cum >= fpca_pve_threshold)[1]
  if (is.na(k_pve)) k_pve <- fp$npc
  K <- max(fpca_min_npc, min(k_pve, fpca_max_npc, fp$npc))
  S <- fp$scores[, 1:K, drop = FALSE]
  achieved <- pve_cum[K]
  cat(sprintf("  FPCA: retained K=%d; achieved cumPVE=%.4f (target=%.4f); FPC1=%.1f%%\n",
              K, achieved, fpca_pve_threshold, 100 * pve_ind[1]))
  
  # --- NEW: per-subject PC1 + FPC1 + FPC2 scores (Supp Fig 19) ---
  fpc1 <- fp$scores[, 1]
  fpc2 <- if (ncol(fp$scores) >= 2) fp$scores[, 2] else rep(NA_real_, nrow(fp$scores))
  metric_scores <- data.frame(
    ID = df_m$ID, metric_key = m, metric_label = metric_labels[m],
    FPC1 = fpc1, FPC2 = fpc2, stringsAsFactors = FALSE
  )
  metric_scores <- merge(metric_scores, pc_scores_df, by = "ID", all.x = TRUE)
  score_rows[[m]] <- metric_scores
  
  # Designs
  X_full       <- model.matrix(reformulate(c(forced_covariates, exposures)), data = df_m)
  X_red_global <- model.matrix(reformulate(forced_covariates), data = df_m)
  
  # Global test
  g <- perm_test_increment(S = S, X_full = X_full, X_red = X_red_global, B = B_perm)
  global_results[[m]] <- tibble(
    metric = m, metric_label = metric_labels[m], K = K,
    achieved_cumPVE = achieved, T_global = g$T_obs,
    p_global_perm = g$p_perm, N = nrow(df_m)
  )
  
  # Incremental tests (each exposure uniquely)
  inc_tbl <- map_dfr(exposures, function(xj) {
    X_red <- model.matrix(reformulate(c(forced_covariates, setdiff(exposures, xj))), data = df_m)
    res <- perm_test_increment(S = S, X_full = X_full, X_red = X_red, B = B_perm)
    tibble(metric = m, metric_label = metric_labels[m],
           exposure = xj, exposure_label = exposure_labels[xj],
           K = K, achieved_cumPVE = achieved, T_inc = res$T_obs, p_inc_perm = res$p_perm)
  }) %>% arrange(p_inc_perm)
  perm_results_all[[m]] <- inc_tbl
  
  # Latent severity axis (exposure PC1) in the same FPCA score space
  exp_scaled_subset <- scale(df_m[, exposures, drop = FALSE],
                             center = attr(exp_for_pca, "scaled:center"),
                             scale  = attr(exp_for_pca, "scaled:scale"))
  sev_pc1 <- as.numeric(exp_scaled_subset %*% pca_exp$rotation[, 1])
  df_sev <- df_m %>% mutate(severity_pc1 = sev_pc1)
  X_full_sev <- model.matrix(reformulate(c(forced_covariates, "severity_pc1")), data = df_sev)
  X_red_sev  <- model.matrix(reformulate(forced_covariates), data = df_sev)
  sev_res <- perm_test_increment(S = S, X_full = X_full_sev, X_red = X_red_sev, B = B_perm)
  latent_results[[m]] <- tibble(
    metric = m, metric_label = metric_labels[m], K = K,
    achieved_cumPVE = achieved, T_severity = sev_res$T_obs,
    p_severity_perm = sev_res$p_perm, N = nrow(df_m),
    exposure_pc1_pve = round(pve[1] * 100, 1)
  )
}

# ==============================================================================
# SAVE COMBINED RESULTS
# ==============================================================================

global_tbl <- bind_rows(global_results) %>%
  mutate(p_global_perm = signif(p_global_perm, 3)) %>% arrange(p_global_perm)
inc_tbl_all <- bind_rows(perm_results_all) %>%
  mutate(p_inc_perm = signif(p_inc_perm, 3)) %>% arrange(metric, p_inc_perm)
latent_tbl <- bind_rows(latent_results) %>%
  mutate(p_severity_perm = signif(p_severity_perm, 3)) %>% arrange(p_severity_perm)

write.csv(global_tbl,  file.path(output_dir, "perm_global_exposure_set_by_metric.csv"),  row.names = FALSE)
write.csv(inc_tbl_all, file.path(output_dir, "perm_incremental_exposure_by_metric.csv"), row.names = FALSE)
write.csv(latent_tbl,  file.path(output_dir, "perm_latent_severity_by_metric.csv"),      row.names = FALSE)

# NEW: per-subject PC1 + FPC scores for Supp Fig 19
if (length(score_rows) > 0) {
  scores_long <- bind_rows(score_rows)
  write.csv(scores_long, file.path(output_dir, "PC1_FPC_scores.csv"), row.names = FALSE)
}

cat("\nSaved:\n")
cat("  perm_global_exposure_set_by_metric.csv\n")
cat("  perm_incremental_exposure_by_metric.csv\n")
cat("  perm_latent_severity_by_metric.csv\n")
cat("  PC1_FPC_scores.csv\n")

cat("\n========================================\n")
cat("Note 7 analysis: DONE\n")
cat(sprintf("Outputs in: %s\n", output_dir))
cat("========================================\n")