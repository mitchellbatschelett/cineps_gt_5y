% ============================================================
% ONE_apply_mu_norm.m
%
% PURPOSE: Apply the SIFT2 proportionality coefficient (mu) to each
%          raw SIFT2-weighted DK connectome, producing mu-weighted
%          matrices suitable for inter-subject comparison of absolute
%          edge weights (Smith et al. 2022).
%
%          For each subject:
%              W_mu_weighted = mu * W_raw
%          where W_raw is the 84x84 raw SIFT2-weighted matrix and mu
%          is the scalar proportionality coefficient from tcksift2.
%
% INPUTS:  For each subject in IN_DIR (subject-subfolder layout):
%            - IN_DIR/<subject>/<subject>_dk.csv
%            - IN_DIR/<subject>/<subject>_sift2_mu.txt
%          (Both produced by FOUR_tractography_and_connectome.sh)
%
% OUTPUTS: For each subject in OUT_DIR (subject-subfolder layout):
%            - OUT_DIR/<subject>/<subject>_mu_weighted.csv
% ============================================================

clear; clc;

% ============================================================
% CONFIG -- edit before running
% ============================================================
% Paths are relative to the repository root; run MATLAB from the repo root
% (or set an absolute path here if launching from elsewhere).
IN_DIR  = 'data/connectomes/raw';         % parent dir containing per-subject subfolders
OUT_DIR = 'data/connectomes/mu_weighted'; % parent output dir; per-subject subfolders created automatically
SUBLIST = '';     % optional: path to text file with subject IDs (one per line). Leave '' to process all sub-* in IN_DIR.
FORCE   = false;  % set true to overwrite existing mu-weighted outputs

% ============================================================
% Logic below -- do not edit
% ============================================================

if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end

% Resolve subject list
if ~isempty(SUBLIST) && exist(SUBLIST, 'file')
    fid = fopen(SUBLIST); subjects = textscan(fid, '%s'); fclose(fid);
    subjects = subjects{1};
else
    d = dir(fullfile(IN_DIR, 'sub-*'));
    subjects = {d([d.isdir]).name}';
end
fprintf('Found %d subjects to process.\n\n', length(subjects));

n_done = 0;
n_skip = 0;
n_fail = 0;

for i = 1:length(subjects)
    sub = subjects{i};
    in_csv  = fullfile(IN_DIR, sub, [sub '_dk.csv']);
    in_mu   = fullfile(IN_DIR, sub, [sub '_sift2_mu.txt']);
    out_dir_sub = fullfile(OUT_DIR, sub);
    out_csv = fullfile(out_dir_sub, [sub '_mu_weighted.csv']);

    % Idempotency check
    if exist(out_csv, 'file') && ~FORCE
        fprintf('[%s] already processed -- skipping (set FORCE=true to re-run)\n', sub);
        n_skip = n_skip + 1;
        continue;
    end

    % Input validation
    if ~exist(in_csv, 'file')
        fprintf('[%s] SKIPPED: %s not found\n', sub, [sub '_dk.csv']);
        n_fail = n_fail + 1;
        continue;
    end
    if ~exist(in_mu, 'file')
        fprintf('[%s] SKIPPED: %s not found\n', sub, [sub '_sift2_mu.txt']);
        n_fail = n_fail + 1;
        continue;
    end

    try
        W  = readmatrix(in_csv);
        mu = readmatrix(in_mu);

        if ~isequal(size(W), [84 84])
            fprintf('[%s] WARNING: connectome is %dx%d (expected 84x84)\n', sub, size(W,1), size(W,2));
        end
        if ~isscalar(mu)
            fprintf('[%s] SKIPPED: sift2_mu.txt is not a scalar (size %dx%d)\n', sub, size(mu,1), size(mu,2));
            n_fail = n_fail + 1;
            continue;
        end

        if ~isfolder(out_dir_sub), mkdir(out_dir_sub); end
        writematrix(W * mu, out_csv);
        fprintf('[%s] done\n', sub);
        n_done = n_done + 1;
    catch ME
        fprintf('[%s] ERROR -- %s\n', sub, ME.message);
        n_fail = n_fail + 1;
    end
end

fprintf('\n============================================================\n');
fprintf('Summary: %d processed, %d skipped, %d failed\n', n_done, n_skip, n_fail);
fprintf('============================================================\n');