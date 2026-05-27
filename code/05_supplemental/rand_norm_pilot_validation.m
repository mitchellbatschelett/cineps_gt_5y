% ============================================================
% rand_norm_pilot_validation.m
%
% PURPOSE: Randomization-normalized GT metrics --- pilot validation with
%          1000 null models on 15 VPT + 15 FT subjects (seed=42), to
%          compare against the 100-null estimates used in the main
%          analysis and assess null-model stability.
%
% INPUTS:  For each selected subject:
%            - IN_DIR/<subject>/<subject>_mu_weighted.csv
%              (same layout as ONE_compute_graph_theory_metrics.m)
%          COHORT_XLSX (ID + Group columns) defines the eligible pool.
%          Random selection uses rng(42, 'twister'): 15 VPT (Group==1)
%          plus 15 FT (Group==0), drawn in cohort-xlsx order.
%
% OUTPUTS: Per subject, written to OUT_DIR:
%            - <subject>_pilot_1000null.mat
%              (fields: subj_id, subj_group, densities,
%                       wei_GE, wei_ACC, wei_L, wei_SW)
%
% USAGE:   Run per-subject from the LSF job array
%          (see submit_rand_norm_pilot.sh):
%            matlab -nodisplay -nosplash -batch ...
%              "rand_norm_pilot_validation(${LSB_JOBINDEX})"
%          where ${LSB_JOBINDEX} is in [1, 30].
%
% DOWNSTREAM: merge_pilot_validation.m combines the 30 .mat files into
%             data/intermediate/pilot_1000null_merged.xlsx, which is
%             consumed by code/05_supplemental/build_supplement.R
%             Section N2 (renders Supplementary Figure 3, Sup Table 2).
%
% REQUIRES: MATLAB R2025b (or similar), Brain Connectivity Toolbox.
% ============================================================

function rand_norm_pilot_validation(subj_idx)

if nargin < 1
    error('Must supply subject index (1-30).');
end

% ============================================================
% CONFIG -- edit before running
% ============================================================
IN_DIR       = 'data/connectomes/mu_weighted';
COHORT_XLSX  = 'data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx';
OUT_DIR      = 'data/intermediate/pilot_validation_results';
BCT_PATH     = '';   % set to BCT install path if not already on the MATLAB path

% --- Density grid (matches main analysis) ---
DENSITY_START = 0.0025;
DENSITY_STEP  = 0.0025;
DENSITY_END   = 1.0;

% --- Null model settings ---
N_RAND        = 1000;   % null networks per subject per threshold (pilot)
WEI_BIN_SWAPS = 10;
WEI_FREQ      = 1;

% --- Charpath settings (match main analysis) ---
CHARPATH_DIAG = 0;
CHARPATH_INF  = 0;

% --- Selection seed (defines which 30 subjects make up the pilot) ---
SELECTION_SEED = 42;
N_VPT_SAMPLE   = 15;
N_FT_SAMPLE    = 15;

% --- Parallelization ---
N_CORES = feature('numcores');

% --- Expected matrix size (sanity check) ---
EXPECTED_N_NODES = 84;

% ============================================================
% Logic below -- do not edit
% ============================================================

% --- BCT path ---
if ~isempty(BCT_PATH)
    addpath(BCT_PATH);
end

% --- Ensure output directory exists ---
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end

% --- Density grid ---
densities = DENSITY_START : DENSITY_STEP : DENSITY_END;
nThresh   = numel(densities);

% --- Read cohort manifest ---
if ~exist(COHORT_XLSX, 'file')
    error('Cohort manifest not found: %s', COHORT_XLSX);
end
manifest = readtable(COHORT_XLSX);
if ~all(ismember({'ID', 'Group'}, manifest.Properties.VariableNames))
    error('Cohort manifest must contain ID and Group columns.');
end

% Eligible IDs by group, in manifest order.
% Handle either numeric or string ID columns (Excel may store IDs as text).
if iscell(manifest.ID)
    ids_num = cellfun(@(x) str2double(x), manifest.ID);
elseif isstring(manifest.ID) || ischar(manifest.ID)
    ids_num = str2double(string(manifest.ID));
else
    ids_num = manifest.ID;
end
if any(isnan(ids_num))
    error('Could not parse all IDs as numeric. Check the ID column in %s.', COHORT_XLSX);
end
vpt_ids = ids_num(manifest.Group == 1);
ft_ids  = ids_num(manifest.Group == 0);

if numel(vpt_ids) < N_VPT_SAMPLE
    error('Manifest has only %d VPT subjects; need %d.', numel(vpt_ids), N_VPT_SAMPLE);
end
if numel(ft_ids) < N_FT_SAMPLE
    error('Manifest has only %d FT subjects; need %d.', numel(ft_ids), N_FT_SAMPLE);
end

% Random selection with fixed seed --- reproducible across runs
rng(SELECTION_SEED, 'twister');
vpt_pick = randperm(numel(vpt_ids), N_VPT_SAMPLE);
ft_pick  = randperm(numel(ft_ids),  N_FT_SAMPLE);

selected_ids    = [vpt_ids(vpt_pick); ft_ids(ft_pick)];
selected_groups = [ones(N_VPT_SAMPLE, 1); zeros(N_FT_SAMPLE, 1)];
n_selected      = numel(selected_ids);   % 30

fprintf('\n=========================================================\n');
fprintf('  PILOT VALIDATION -- 1000 NULL MODELS\n');
fprintf('=========================================================\n');
fprintf('  Selected:        %d VPT + %d FT = %d subjects (seed=%d)\n', ...
        N_VPT_SAMPLE, N_FT_SAMPLE, n_selected, SELECTION_SEED);
fprintf('  Subject index:   %d / %d\n', subj_idx, n_selected);
fprintf('  Null models:     %d per threshold\n', N_RAND);
fprintf('  charpath_inf:    %d (matches main analysis)\n', CHARPATH_INF);
fprintf('  Cores:           %d\n', N_CORES);
fprintf('=========================================================\n\n');

% --- Validate subject index ---
if subj_idx < 1 || subj_idx > n_selected
    error('subj_idx %d is out of range [1, %d]', subj_idx, n_selected);
end

% --- Pick this subject ---
subj_id_num = selected_ids(subj_idx);
subj_group  = selected_groups(subj_idx);
subj_id     = sprintf('sub-%d', subj_id_num);

filepath = fullfile(IN_DIR, subj_id, [subj_id '_mu_weighted.csv']);

fprintf('Processing subject %d/%d: %s  (group=%d)\n', ...
        subj_idx, n_selected, subj_id, subj_group);
fprintf('File: %s\n\n', filepath);

if ~isfile(filepath)
    error('Connectome file not found: %s', filepath);
end

% --- Idempotency check ---
out_file = fullfile(OUT_DIR, sprintf('%s_pilot_1000null.mat', subj_id));
if isfile(out_file)
    fprintf('Output already exists: %s\n', out_file);
    fprintf('Skipping --- delete file and resubmit to rerun.\n');
    return;
end

% --- Load and validate matrix ---
W_raw = readmatrix(filepath);
if ~isequal(size(W_raw), [EXPECTED_N_NODES, EXPECTED_N_NODES])
    error('Dimension mismatch for %s: got %dx%d, expected %dx%d', ...
          subj_id, size(W_raw,1), size(W_raw,2), ...
          EXPECTED_N_NODES, EXPECTED_N_NODES);
end

W_raw = (W_raw + W_raw.') / 2;
n = size(W_raw, 1);
W_raw(1:n+1:end) = 0;
W_raw(W_raw < 0) = 0;
fprintf('Matrix loaded and validated: %dx%d\n\n', n, n);

% --- Parallel pool ---
if isempty(gcp('nocreate'))
    parpool('local', N_CORES);
end
fprintf('Parallel pool: %d workers\n\n', gcp().NumWorkers);

% --- Preallocate ---
wei_GE  = NaN(1, nThresh);
wei_ACC = NaN(1, nThresh);
wei_L   = NaN(1, nThresh);
wei_SW  = NaN(1, nThresh);

% --- Density loop ---
subj_tic = tic;

for t = 1:nThresh
    p = densities(t);

    W_thr = threshold_proportional(W_raw, p);
    if ~any(W_thr(:)), continue; end

    % Observed weighted metrics
    obs_wei_GE  = efficiency_wei(W_thr);
    Ci_wei      = clustering_coef_wu(W_thr);
    obs_wei_ACC = mean(Ci_wei(~isinf(Ci_wei) & ~isnan(Ci_wei)));
    L_mat_wei   = weight_conversion(W_thr, 'lengths');
    D_wei       = distance_wei(L_mat_wei);
    obs_wei_L   = charpath(D_wei, CHARPATH_DIAG, CHARPATH_INF);

    % Null distributions
    null_wei_GE  = NaN(N_RAND, 1);
    null_wei_ACC = NaN(N_RAND, 1);
    null_wei_L   = NaN(N_RAND, 1);

    W_loc = W_thr;
    p_bin = WEI_BIN_SWAPS;
    p_frq = WEI_FREQ;
    p_cd  = CHARPATH_DIAG;
    p_ci  = CHARPATH_INF;

    parfor r = 1:N_RAND
        W_null          = null_model_und_sign(W_loc, p_bin, p_frq);
        null_wei_GE(r)  = efficiency_wei(W_null);
        Ci_nw           = clustering_coef_wu(W_null);
        null_wei_ACC(r) = mean(Ci_nw(~isinf(Ci_nw) & ~isnan(Ci_nw)));
        Ln              = weight_conversion(W_null, 'lengths');
        Dn              = distance_wei(Ln);
        null_wei_L(r)   = charpath(Dn, p_cd, p_ci);
    end

    mn_GE  = mean(null_wei_GE,  'omitnan');
    mn_ACC = mean(null_wei_ACC, 'omitnan');
    mn_L   = mean(null_wei_L,   'omitnan');

    if mn_GE  ~= 0, wei_GE(t)  = obs_wei_GE  / mn_GE;  end
    if mn_ACC ~= 0, wei_ACC(t) = obs_wei_ACC / mn_ACC; end
    if mn_L   ~= 0, wei_L(t)   = obs_wei_L   / mn_L;   end
    if mn_ACC ~= 0 && mn_L ~= 0
        gamma  = obs_wei_ACC / mn_ACC;
        lambda = obs_wei_L   / mn_L;
        if lambda ~= 0, wei_SW(t) = gamma / lambda; end
    end

    if mod(t, 50) == 0 || t == nThresh
        fprintf('  Threshold %3d/%d (density=%.4f) --- elapsed: %.1f min\n', ...
                t, nThresh, p, toc(subj_tic)/60);
    end
end

% --- Cleanup and save ---
pool = gcp('nocreate');
if ~isempty(pool), delete(pool); end

save(out_file, ...
     'subj_id', 'subj_group', 'densities', ...
     'wei_GE', 'wei_ACC', 'wei_L', 'wei_SW');

fprintf('\nSaved: %s\n', out_file);
fprintf('Total time: %.1f minutes\n', toc(subj_tic)/60);
fprintf('=========================================================\n\n');

end