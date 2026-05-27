% ============================================================
% merge_pilot_validation.m
%
% PURPOSE: Combine per-subject .mat outputs from rand_norm_pilot_validation.m
%          into a single Excel file matching the column-name convention of
%          the main analytic xlsx (data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx).
%
% INPUTS:  IN_DIR/sub-<id>_pilot_1000null.mat (one per pilot subject)
%
% OUTPUTS: OUT_FILE with columns:
%            - ID (numeric, no sub- prefix)
%            - Group (1 = VPT, 0 = FT)
%            - rand_norm_wei_GE_<dens>   (400 cols)
%            - rand_norm_wei_ACC_<dens>  (400 cols)
%            - rand_norm_wei_L_<dens>    (400 cols)
%            - rand_norm_wei_SW_<dens>   (400 cols)
%
% USAGE:   Run after all 30 pilot LSF array jobs complete:
%            matlab -nodisplay -nosplash -batch "merge_pilot_validation"
%
% CONSUMED BY: code/05_supplemental/build_supplement.R, Section N2
%              (renders Supplementary Figure 3 and Sup Table 2).
%
% REQUIRES: MATLAB.
% ============================================================

clear; clc;

% ============================================================
% CONFIG -- edit before running
% ============================================================
IN_DIR   = 'data/intermediate/pilot_validation_results';
OUT_FILE = 'data/intermediate/pilot_1000null_merged.xlsx';
FORCE    = false;   % set true to overwrite existing OUT_FILE

% ============================================================
% Logic below -- do not edit
% ============================================================

% Idempotency
if exist(OUT_FILE, 'file') && ~FORCE
    fprintf('Output already exists: %s\n', OUT_FILE);
    fprintf('Set FORCE=true to overwrite. Exiting.\n');
    return;
end

% Ensure output directory exists
[out_dir, ~, ~] = fileparts(OUT_FILE);
if ~isfolder(out_dir), mkdir(out_dir); end

% --- Locate per-subject .mat files ---
mat_files = dir(fullfile(IN_DIR, 'sub-*_pilot_1000null.mat'));
n_files   = numel(mat_files);

if n_files == 0
    error('No pilot .mat files found in: %s', IN_DIR);
end

fprintf('Found %d pilot result files.\n\n', n_files);

% --- Load all subjects ---
rows = cell(n_files, 1);
for i = 1:n_files
    rows{i} = load(fullfile(IN_DIR, mat_files(i).name));
end

densities = rows{1}.densities;
nThresh   = numel(densities);

% --- Build density column labels (zero-padded 6-char, e.g. 000.25, 100.00) ---
density_labels = cell(1, nThresh);
for t = 1:nThresh
    pct = densities(t) * 100;
    if pct >= 100
        density_labels{t} = '100.00';
    else
        density_labels{t} = sprintf('%06.2f', pct);
    end
end

% --- ID and Group columns ---
% IDs are stored as 'sub-<id>'; strip the sub- prefix and cast to numeric
% to match the convention in cohort_171VPT_45FT_postVQC.xlsx.
IDs    = cellfun(@(r) str2double(strrep(r.subj_id, 'sub-', '')), rows);
Groups = cellfun(@(r) r.subj_group, rows);

T = table(IDs, Groups, 'VariableNames', {'ID', 'Group'});

% --- Density-resolved columns ---
for t = 1:nThresh
    dl = density_labels{t};
    T.(sprintf('rand_norm_wei_GE_%s',  dl)) = cellfun(@(r) r.wei_GE(t),  rows);
    T.(sprintf('rand_norm_wei_ACC_%s', dl)) = cellfun(@(r) r.wei_ACC(t), rows);
    T.(sprintf('rand_norm_wei_L_%s',   dl)) = cellfun(@(r) r.wei_L(t),   rows);
    T.(sprintf('rand_norm_wei_SW_%s',  dl)) = cellfun(@(r) r.wei_SW(t),  rows);
end

% --- Write ---
writetable(T, OUT_FILE);

fprintf('=========================================================\n');
fprintf('  MERGE COMPLETE\n');
fprintf('=========================================================\n');
fprintf('  Subjects:    %d\n', n_files);
fprintf('  Thresholds:  %d\n', nThresh);
fprintf('  Columns:     %d\n', width(T));
fprintf('  Output:      %s\n', OUT_FILE);
fprintf('=========================================================\n');