#!/bin/bash
# ============================================================
# submit_rand_norm_pilot.sh
#
# Submits a 30-element LSF job array to run rand_norm_pilot_validation.m
# on each of the 30 pilot subjects. Each task processes one subject with
# 1000 null models per density threshold.
#
# Usage (from the repo root):
#   bsub < code/05_supplemental/submit_rand_norm.sh
#
# Requires:
#   - MATLAB module on the HPC
#   - Brain Connectivity Toolbox available on the MATLAB path (or set
#     BCT_PATH in rand_norm_pilot_validation.m before running)
#   - Repo cloned at $REPO_ROOT below
# ============================================================

#BSUB -J pilot_val[1-30]%1
#BSUB -n 48
#BSUB -M 64000
#BSUB -R "span[hosts=1]"
#BSUB -W 24:00
#BSUB -q normal
#BSUB -o logs/pilot_val_%I.out
#BSUB -e logs/pilot_val_%I.err

# --- Repo root (edit if running from a different location) ---
REPO_ROOT="${REPO_ROOT:-$(pwd)}"   # defaults to current dir; run from repo root, or set REPO_ROOT

cd "$REPO_ROOT" || { echo "Repo root not found: $REPO_ROOT" >&2; exit 1; }

mkdir -p logs
mkdir -p data/intermediate/pilot_validation_results

module purge
module load matlab/2025b

matlab -nodisplay -nosplash -batch \
    "addpath('code/05_supplemental'); rand_norm_pilot_validation(${LSB_JOBINDEX})"