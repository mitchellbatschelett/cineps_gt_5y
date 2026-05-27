% ============================================================
% network_fragmentation.m
%
% PURPOSE: Compute the percentage of subjects with a single connected
%          component as a function of network density (0.25%-100% in
%          0.25% steps), separately for All / VPT / FT groups.
%
%          The minimum density at which all subjects exhibit a single
%          connected component (here: 11%) establishes the empirical
%          lower bound for the FDA density sweep used throughout the
%          manuscript.
%
% INPUTS:  For each subject listed in COHORT_XLSX (column: ID, Group):
%            - IN_DIR/<subject>/<subject>_mu_weighted.csv
%              (output of ONE_apply_mu_norm.m, same layout as
%              ONE_compute_graph_theory_metrics.m)
%          COHORT_XLSX is read for ID (numeric) and Group (1=VPT, 0=FT)
%          to label each subject's group.
%
% OUTPUTS: One CSV at OUT_FILE with columns:
%            - density_pct
%            - pct_single_comp_all
%            - pct_single_comp_vpt
%            - pct_single_comp_ft
%
% CONSUMED BY: code/05_supplemental/build_supplement.R, Section N3
%              (renders Supplementary Figure 4).
%
% REQUIRES: MATLAB (no BCT functions needed; BFS implemented inline).
% ============================================================

clear; clc; close all;

% ============================================================
% CONFIG -- edit before running
% ============================================================
IN_DIR       = 'data/connectomes/mu_weighted';                          % parent dir with per-subject subfolders
COHORT_XLSX  = 'data/analysis_ready/cohort_171VPT_45FT_postVQC.xlsx';   % manifest: ID + Group columns
OUT_FILE     = 'data/intermediate/fragmentation_by_density.csv';
FORCE        = false;   % set true to overwrite existing OUT_FILE

% --- Density grid (matches FDA analysis: 0.25% steps from 0.25 to 100) ---
DENSITY_START = 0.25;
DENSITY_STEP  = 0.25;
DENSITY_END   = 100;

% --- Expected matrix size (sanity check) ---
EXPECTED_N_NODES = 84;

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

% Density grid
densities   = DENSITY_START : DENSITY_STEP : DENSITY_END;
n_densities = numel(densities);

% --- Read cohort manifest ---
if ~exist(COHORT_XLSX, 'file')
    error('Cohort manifest not found: %s', COHORT_XLSX);
end
manifest = readtable(COHORT_XLSX);
if ~all(ismember({'ID', 'Group'}, manifest.Properties.VariableNames))
    error('Cohort manifest must contain ID and Group columns.');
end

% Build subject IDs in BIDS form (sub-<id>) and group vector.
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
subject_ids = arrayfun(@(x) sprintf('sub-%d', x), ids_num, 'UniformOutput', false);
group_label = manifest.Group;   % 1 = VPT, 0 = FT
n_subjects  = numel(subject_ids);
n_vpt = sum(group_label == 1);
n_ft  = sum(group_label == 0);

fprintf('\n=========================================================\n');
fprintf('  NETWORK FRAGMENTATION ANALYSIS\n');
fprintf('=========================================================\n');
fprintf('  Cohort:       %d subjects (%d VPT, %d FT)\n', n_subjects, n_vpt, n_ft);
fprintf('  Connectomes:  %s\n', IN_DIR);
fprintf('  Densities:    %d thresholds (%.2f : %.2f : %.2f)\n', ...
        n_densities, DENSITY_START, DENSITY_STEP, DENSITY_END);
fprintf('  Output:       %s\n', OUT_FILE);
fprintf('=========================================================\n\n');

% --- Load all connectivity matrices ---
fprintf('=== Loading Connectivity Matrices ===\n');
n_regions    = EXPECTED_N_NODES;
all_matrices = zeros(n_regions, n_regions, n_subjects);
valid_mask   = true(n_subjects, 1);

for s = 1:n_subjects
    sub    = subject_ids{s};
    in_csv = fullfile(IN_DIR, sub, [sub '_mu_weighted.csv']);

    if ~exist(in_csv, 'file')
        fprintf('  [missing] %s\n', in_csv);
        valid_mask(s) = false;
        continue;
    end

    mat = readmatrix(in_csv);
    if ~isequal(size(mat), [n_regions, n_regions])
        fprintf('  [skip] %s: dim %dx%d, expected %dx%d\n', sub, ...
                size(mat,1), size(mat,2), n_regions, n_regions);
        valid_mask(s) = false;
        continue;
    end

    % Symmetrize, zero diagonal
    mat = (mat + mat') / 2;
    mat(logical(eye(n_regions))) = 0;
    all_matrices(:, :, s) = mat;
end

% Trim invalid subjects
subject_ids  = subject_ids(valid_mask);
group_label  = group_label(valid_mask);
all_matrices = all_matrices(:, :, valid_mask);
n_subjects   = numel(subject_ids);
fprintf('  Loaded: %d/%d subjects\n', n_subjects, numel(valid_mask));

% --- Extract upper-triangle edge weights per subject ---
n_edges  = n_regions * (n_regions - 1) / 2;
tri_idx  = find(triu(ones(n_regions), 1));
edge_weights_all = zeros(n_edges, n_subjects);
for s = 1:n_subjects
    mat = all_matrices(:, :, s);
    edge_weights_all(:, s) = mat(tri_idx);
end
[edge_row, edge_col] = ind2sub([n_regions, n_regions], tri_idx);

% --- Main loop: per-subject, per-density connectedness ---
fprintf('\n=== Computing Connectedness Across %d Densities ===\n', n_densities);

% is_single_component(s, d) = true if at density d subject s's network has
% a single connected component among all active (degree > 0) nodes.
is_single_component = false(n_subjects, n_densities);

for s = 1:n_subjects
    subj_weights = edge_weights_all(:, s);
    [sorted_weights, sort_order] = sort(subj_weights, 'descend');

    % Build adjacency incrementally as we sweep density upward
    adj      = zeros(n_regions, n_regions);
    edge_ptr = 0;

    for d = 1:n_densities
        n_keep = round(densities(d) / 100 * n_edges);

        % Add edges from edge_ptr+1 up to n_keep
        while edge_ptr < n_keep && edge_ptr < n_edges
            edge_ptr = edge_ptr + 1;

            % Stop at the first zero-weight edge (sorted descending)
            if sorted_weights(edge_ptr) == 0
                edge_ptr = n_keep;
                break;
            end

            ei = sort_order(edge_ptr);
            r  = edge_row(ei);
            c  = edge_col(ei);
            adj(r, c) = 1;
            adj(c, r) = 1;
        end

        % Active nodes: degree > 0
        node_degree  = sum(adj, 2);
        active_nodes = node_degree > 0;
        n_active     = sum(active_nodes);

        % Edge case: zero active nodes -> vacuously single component
        if n_active == 0
            is_single_component(s, d) = true;
            continue;
        end

        % BFS to count components among active nodes
        visited = false(n_regions, 1);
        n_comp  = 0;

        for start_node = 1:n_regions
            if ~active_nodes(start_node) || visited(start_node)
                continue;
            end

            n_comp = n_comp + 1;
            queue  = start_node;
            visited(start_node) = true;

            while ~isempty(queue)
                node     = queue(1);
                queue(1) = [];
                neighbors = find(adj(node, :) > 0);
                for ni = 1:length(neighbors)
                    nb = neighbors(ni);
                    if ~visited(nb) && active_nodes(nb)
                        visited(nb) = true;
                        queue(end + 1) = nb;  %#ok<AGROW>
                    end
                end
            end
        end

        is_single_component(s, d) = (n_comp == 1);
    end

    if mod(s, 25) == 0
        fprintf('  %d/%d subjects done\n', s, n_subjects);
    end
end
fprintf('  %d/%d subjects done\n', n_subjects, n_subjects);

% --- Group-level percentages per density ---
vpt_mask = group_label == 1;
ft_mask  = group_label == 0;

pct_single_comp_all = mean(is_single_component(:,        :), 1) * 100;
pct_single_comp_vpt = mean(is_single_component(vpt_mask, :), 1) * 100;
pct_single_comp_ft  = mean(is_single_component(ft_mask,  :), 1) * 100;

% Identify the all-connected threshold (matches Sup Fig 4 vertical dashed line)
all_connected_idx = find(pct_single_comp_all >= 100, 1, 'first');
if ~isempty(all_connected_idx)
    fprintf('\n  All subjects single-component at density: %.2f%%\n', ...
            densities(all_connected_idx));
end

% --- Save CSV ---
density_table = table( ...
    densities', pct_single_comp_all', pct_single_comp_vpt', pct_single_comp_ft', ...
    'VariableNames', {'density_pct', 'pct_single_comp_all', ...
                      'pct_single_comp_vpt', 'pct_single_comp_ft'});

writetable(density_table, OUT_FILE);

fprintf('\n=========================================================\n');
fprintf('  COMPLETE\n');
fprintf('=========================================================\n');
fprintf('  Output: %s\n', OUT_FILE);
fprintf('=========================================================\n');