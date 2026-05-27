# Reproducibility repository for:

> **Title:** The effects of very premature birth and prematurity-related morbidities on brain structural connectivity at early school age
>
> **Authors:** Mitchell A. Batschelett, Armin Allahverdy, Brady J. Williamson, Mekibib Altaye, Beth Kline-Fath, Jean Tkach, Weihong Yuan, & Nehal A. Parikh

This repository contains the code, connectomes, derived analytic tables, and saved statistical results needed to reproduce every table and figure in the manuscript and supplement.

---

## Two ways to use this repository

**Path A — regenerate tables and figures from saved results (minutes).**
The repository contains the saved statistical outputs in `results/`, so you can rebuild every main and supplementary table and figure without re-running any of the computationally expensive FDA modeling. This is the recommended path for most users and for reproducibility review. Two scripts do all of it:

```bash
Rscript code/04_results_figures_tables/build_main_results.R   # Tables 1-3, Figures 2-5
Rscript code/05_supplemental/build_supplement.R               # Supp Tables 1-7, Supp Figures 1-20
```

Both scripts locate the repository root automatically (see "How paths work" below) and write into `tables/` and `figures/`.

**Path B — re-run the full statistical pipeline from the connectomes (days to a week).**
If you want to regenerate the saved results themselves, you can run the analysis pipeline starting from the shared connectomes. This is computationally expensive. Depending on computational resources, the FDA group-difference models take hours each, the small-worldness (scat-family) models take up to a week, and stability selection takes roughly a day per metric. We ran these jobs on an HPC cluster; the `submit_*.lsf` / `submit_*.sh` scripts are provided as illustrative cluster launchers and will need site-specific edits.

You do **not** need Path B to reproduce the manuscript's tables and figures.

### Running on a cluster (Path B detail)

The functional data analyses were run as LSF array jobs. Each launcher (`code/03_statistical_analysis/submit_FDA_array_*.lsf`) sets up the environment, changes into the statistical analysis directory, and calls the corresponding *runs* driver with an array index, e.g.:

```bash
SCRIPT_DIR=/path/to/code/03_statistical_analysis   # <-- edit to your path
cd $SCRIPT_DIR
Rscript TWO_fda_group_difference_runs.R $LSB_JOBINDEX
```

Before submitting, edit each launcher for your site:

- set `SCRIPT_DIR` to your unique `code/03_statistical_analysis/` path;
- set the `#BSUB -e` / `#BSUB -o` log paths to writable locations;
- adjust the scheduler directives (queue, wall time, cores, memory) and the `module load R/...` line to your cluster.

Because the launcher `cd`s into `SCRIPT_DIR`, the runs scripts find the repository root automatically from there — no path inside the R scripts needs editing. The same pattern applies to `code/05_supplemental/submit_rand_norm.sh`, which honors a `REPO_ROOT` environment variable (defaulting to the current directory).

Note that the run tables and per-job cost are documented in the header of each launcher and each runs driver. In particular, the small-worldness (scat-family) group-difference runs are far slower than the others and should be submitted with extended wall time.

If a cluster is not available or you wish to run the analyses locally, the pipeline scripts themselves can be run directly. Configurations for each run can be tuned in specific code blocks near the top of the scripts.

---

## Included data

**Included:** the connectomes (84 × 84 Desikan-Killiany matrices), the derived analytic cohort tables, all saved statistical results, and all code from connectome construction onward.

**Not included:** no raw or preprocessed MRI, and no upstream image-processing code. The anatomical preprocessing, motion QC, diffusion preprocessing, and tractography stages are described in the manuscript Methods but their code and data are not included with this repository. Analysis here begins from the connectomes. This is due to the ongoing nature of the study, as well as privacy agreements with participants and their families.

### Shared connectomes

The raw connectomes in `data/connectomes/raw/<subject>/` (`<subject>_dk.csv` plus `<subject>_sift2_mu.txt`) are the output of the tractography pipeline described in manuscript Methods. The first step in this repository is mu-weighting (`code/01_connectome_construction/`).

---

## Pipeline overview

The repository is organized as a four-stage chain. Within `code/`, directories are numbered in execution order, and scripts within each stage are prefixed `ONE_`, `TWO_`, etc.

| Stage | Directory | What it does | Runs on |
|-------|-----------|--------------|---------|
| 1 | `code/01_connectome_construction/` | Apply the SIFT2 proportionality coefficient (mu) to each raw matrix to generate mu-weighted connectomes (the mu-weighted outputs are also pre-computed and present within `data/connectomes/mu_weighted/`, so this stage is optional) | Local (MATLAB) |
| 2 | `code/02_graph_theory/` | Compute graph metrics across 400 density thresholds + null normalization (MATLAB), then merge with the clinical/demographic table into the analytic cohort files (R) | Local (MATLAB + R) |
| 3 | `code/03_statistical_analysis/` | Density group difference (R); FDA group differences, FDA univariate associations, and FDA stability selection (R, via per-run driver scripts) | Cluster (expensive) |
| 4 | `code/04_results_figures_tables/` and `code/05_supplemental/` | Regenerate every main and supplementary table and figure from saved results | Local (R) |

Within stage 3, each analysis has a **pipeline** script (the analysis itself) and a **runs** driver (selects one run by command-line index, sets parameters, and sources the pipeline). The drivers are what the cluster array jobs call.

---

## How paths work

All R scripts locate the repository root automatically. There should be no paths to edit before running. Each script searches, in order:

1. a `repo_root` you have already set in the R global environment (if any);
2. the script's own file location (when `source()`-d), walking upward;
3. the current working directory, walking upward;

stopping at the first folder that contains both `data/analysis_ready/` and `code/`. If none is found, the script stops with a clear error rather than guessing.

**Local use.** Running either builder with `Rscript` from the repository root, or `source()`-ing it from an R session, both resolve automatically:

```bash
# from the repository root
Rscript code/04_results_figures_tables/build_main_results.R
```

```r
# or from an R session
source("code/04_results_figures_tables/build_main_results.R")

# only needed if auto-detection fails on an unusual setup:
repo_root <- "/full/path/to/repository"
source("code/04_results_figures_tables/build_main_results.R")
```

**Cluster use.** The provided LSF launchers `cd` into `code/03_statistical_analysis/` before calling `Rscript`, so detection succeeds via the working-directory search (step 3 above) — it walks up from the code directory to the repository root. The only thing you must edit in the launchers is the `SCRIPT_DIR` line (and the log paths); see "Running on a cluster" above. If you launch the runs scripts some other way, either `cd` into the repository first or set `repo_root` at the top of the script.

The MATLAB scripts use paths relative to the repository root; run MATLAB from the repository root (or edit the `IN_DIR` / `OUT_FILE` lines at the top of each script).

---

## Software requirements

**R** (developed under R ≥ 4.4.0) with the following packages:

```
refund, mgcv, fda, fda.usc, grpreg,
tidyverse (dplyr, tidyr, tibble, ggplot2), readxl, writexl, openxlsx,
car, cowplot, patchwork, gridExtra, viridis, irr, parallel
```

**MATLAB** (developed under R2025b) with the **Brain Connectivity Toolbox** for the graph-theory metric computation in stage 2. Set `BCT_PATH` at the top of the relevant `.m` scripts if BCT is not already on your MATLAB path. The Parallel Computing Toolbox is also used for graph metric generation and normalization.

---

## Directory map

```
.
├── code/
│   ├── 01_connectome_construction/   mu-weighting (MATLAB)
│   ├── 02_graph_theory/              metric computation (MATLAB) + cohort merge (R)
│   ├── 03_statistical_analysis/      density, FDA group/univariate/stability (R) + cluster launchers
│   ├── 04_results_figures_tables/    build_main_results.R
│   └── 05_supplemental/              build_supplement.R + supplemental analyses
├── data/
│   ├── connectomes/
│   │   ├── raw/<subject>/            <subject>_dk.csv + <subject>_sift2_mu.txt
│   │   └── mu_weighted/<subject>/    <subject>_mu_weighted.csv  (output of stage 1; pre-computed and provided)
│   ├── intermediate/                 per-subject metrics workbook, pilot null validation, fragmentation
│   ├── demographic_clinical/         demo_clinical_only.xlsx  (used by Table 1)
│   └── analysis_ready/               merged analytic cohort files (see Data conventions)
├── results/                          saved statistical outputs (inputs to the builders)
├── figures/                          main/ and supplement/  (regenerated by the builders)
└── tables/                           main/ and supplement/  (regenerated by the builders)
```

---

## Code: what each script does

### `code/01_connectome_construction/`

- **`ONE_apply_mu_norm.m`** — applies the SIFT2 proportionality coefficient (mu) from `<subject>_sift2_mu.txt` to each raw connectome, producing the mu-weighted matrix used by all downstream analyses. The mu-weighted outputs are pre-computed and contained in `data/connectomes/mu_weighted/`, so re-running this script is optional.

### `code/02_graph_theory/`

- **`ONE_compute_graph_theory_metrics.m`** — computes weighted graph theory metrics (global efficiency, average clustering coefficient, mean nodal strength) across 400 proportional density thresholds, plus null-normalized versions of GE, ACC, characteristic path length, and small-worldness. Null models use 100 random networks per subject per threshold (`null_model_und_sign` with `bin_swaps = 10`, `wei_freq = 1`). Output is the per-subject workbook `data/intermediate/graph_theory_metrics_per_subject.xlsx`. This step is the bulk of stage 2's runtime; the workbook is contained, so re-running it is optional.
- **`TWO_merge_clinical_and_metrics.R`** — merges the per-subject metrics workbook with `demo_clinical_only.xlsx` and writes the two analytic cohort files in `data/analysis_ready/`.

### `code/03_statistical_analysis/`

- **`ONE_density_group_difference.R`** — tests group differences in unthresholded connectome density (Welch's *t* + ANCOVA), with a 2-SD outlier sensitivity branch. Produces Figure 1, the density column of Table 2, and Supplementary Table 3.
- **`TWO_fda_pipeline.R` / `TWO_fda_group_difference_runs.R`** — fits density-resolved penalized function-on-scalar regression (PFFR) models for VPT-vs-FT group differences. The pipeline does one model; the runs driver defines the 36 (metric × covariate-configuration) cells and selects one by `$LSB_JOBINDEX`. Uses the `scat` family for the small-worldness metric and gaussian otherwise.
- **`THREE_fda_univariate_pipeline.R` / `THREE_fda_univariate_runs.R`** — fits per-exposure PFFR models within the VPT cohort, fully adjusted for forced covariates (corrected age at MRI, eTIV, sex, relative motion, social risk score). The runs driver defines the 28 (exposure × metric) cells. Used for Supplementary Figures 11-14 and Supplementary Table 6.
- **`FOUR_fda_stability_selection_pipeline.R` / `FOUR_fda_stability_selection_runs.R`** — runs the density-resolved stability-selection procedure for each metric: 100 random 60/40 splits with group LASSO, then re-fits the final PFFR with selected exposures and bootstrap 95% CIs. The runs driver defines the four main-text metric runs plus two ACC sensitivity branches (`high_gba_rem` and `gba_binary`).
- **`submit_FDA_array_group_diff.lsf` / `submit_FDA_array_univariate.lsf`** — illustrative LSF array launchers for the 36 group-difference cells and the 28 univariate cells, respectively. Each launcher documents its run table and per-job resource requests in its header.

### `code/04_results_figures_tables/`

- **`build_main_results.R`** — generates the main-text figures (2-5) and tables (1-3) from `results/` and the analytic cohort files. Reads the saved FDA outputs (per-run summaries plus RData) and writes CSVs to `tables/main/` and PNG+PDF figures to `figures/main/`. This is the script that will reproduce all of the main results presented in the manuscript.

### `code/05_supplemental/`

- **`build_supplement.R`** — generates every Supplementary Figure (1-20) and Supplementary Table (1-7) from saved results and analytic data, using the same path, helper, and naming conventions as `build_main_results.R`. Organized into seven Notes; each Note prints its own progress banner.
- **`note7_shared_severity_analysis.R`** — standalone analysis generator for Supplementary Note 7 (the shared-severity / exposure-PCA framework). Runs the per-metric PFFR residualization → FPCA → Freedman-Lane permutation procedure (B = 1000), writes outputs to `results/post_lasso_shared_severity/`, and exports the per-subject PC1 + FPC1 + FPC2 scores (`PC1_FPC_scores.csv`) so the supplement builder can render Supplementary Figure 19 without re-running the FPCA.
- **`rand_norm_pilot_validation.m`** — pilot validation comparing graph-metric normalization with 100 versus 1000 random networks, on a 30-subject subset. Outputs per-subject `.mat` files to `data/intermediate/pilot_validation_results/`. The number of random networks is one of the key parameters governing null-model precision; this pilot established that 100 is sufficient.
- **`merge_pilot_validation.m`** — combines the 30 per-subject `.mat` outputs from `rand_norm_pilot_validation.m` into the merged workbook `data/intermediate/pilot_1000null_merged.xlsx`, which the supplement builder reads to produce Supplementary Note 2 (null-model validation: Supplementary Table 2, Supplementary Figure 3).
- **`network_fragmentation.m`** — computes the percentage of subjects with a single connected component at each of the 400 density thresholds, used to establish the 11% lower density floor. Output is `data/intermediate/fragmentation_by_density.csv`, which feeds Supplementary Figure 4.
- **`submit_rand_norm.sh`** — illustrative LSF launcher for `rand_norm_pilot_validation.m`. Submits a 30-element array, one job per subject in the pilot subset. Similar usage to the .lsf files used to call the FDA univariate and group difference models.

---

## Data conventions

**Cohort files** (`data/analysis_ready/`): two merged tables produced by `code/02_graph_theory/TWO_merge_clinical_and_metrics.R`.
- `cohort_171VPT_45FT_postVQC.xlsx` — full 216-subject cohort (2818 columns)
- `cohort_171VPT_postVQC.xlsx` — VPT-only, 171 subjects (2817 columns; `Group` dropped)

**Per-density metric columns.** Each cohort file carries 7 metric types (`GE`, `ACC`, `str`, `rand_norm_wei_GE`, `rand_norm_wei_ACC`, `rand_norm_wei_L`, `rand_norm_wei_SW`) across 400 density thresholds, plus `ID` and `den_100.00`. Density labels are zero-padded 5-character strings (`00.25` … `100.00`) for all metrics **except** `str_`, which is unpadded (`str_0.25` … `str_100.00`). Analyses are restricted to the 11–100% density range (the 11% floor is the smallest density at which all subjects form a single connected component).

**Intermediate data** (`data/intermediate/`):

- **`graph_theory_metrics_per_subject.xlsx`** — the per-subject metrics workbook written by `code/02_graph_theory/ONE_compute_graph_theory_metrics.m`. Provided here; re-runnable, or skippable if you use the already-merged cohort files.
- **`pilot_1000null_merged.xlsx`** — merged 30-subject pilot output from `rand_norm_pilot_validation.m` + `merge_pilot_validation.m`. Used by Supplementary Note 2 for null-model convergence checks.
- **`pilot_validation_results/*.mat`** — the 30 per-subject pilot `.mat` files that feed into the merged workbook.
- **`fragmentation_by_density.csv`** — single-connected-component percentage at each of the 400 density thresholds (from `network_fragmentation.m`); used to establish the 11% floor.

**Saved statistical results** (`results/`) — every input the two builders read.

- `density_group_difference/density_sensitivity_results.xlsx` — outputs from `ONE_density_group_difference.R`.
- `fda_group_differences/<metric>_FDA_Group_<covariates>_11-100_fullsample/` — one directory per cell of the group-difference grid, each containing the fitted `*_FDA_results.RData`, PFFR diagnostics and residuals PNGs, and per-predictor + summary CSVs.
- `fda_univariate/<metric>/<metric>_univariate_<exposure>_11-100/` — one directory per (metric × exposure) cell of the fully-adjusted univariate grid.
- `fda_stability_selection/<metric>_stabsel_11-100/` — one directory per main-text metric run, plus `ACC_sensitivity/` for the two ACC sensitivity branches. Each contains the `*_stabsel_results.RData`, selection matrix (CSV + RDS), per-exposure final-fit beta(d) PNGs, and the stability-selection frequency plot.
- `post_lasso_shared_severity/` — outputs from `note7_shared_severity_analysis.R`: exposure clusters, correlation matrix, PCA loadings/variance, VIFs, the per-subject `PC1_FPC_scores.csv`, and the permutation results (global, incremental, and latent-severity).

**Figures and tables** (`figures/`, `tables/`) — both populated by the builders. `figures/main/` holds Figures 1-5 (PNG + PDF), `figures/supplement/` holds Supplementary Figures 1-20 (PNG + PDF), `tables/main/` holds Tables 1-3 (CSV), and `tables/supplement/` holds Supplementary Tables 1-7 (CSV).

---