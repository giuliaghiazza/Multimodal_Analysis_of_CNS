% GROUP STATISTICS — Task vs. Rest: GFP, HFP, AI
% Level metrics plus temporal-stability metrics (10s non-overlapping windows).

clear; clc;

base_dir    = fileparts(mfilename('fullpath'));
results_dir = fullfile(base_dir, 'RESULTS');

load(fullfile(results_dir, 'Group_GFP_HFP_Results.mat'));
load(fullfile(results_dir, 'Group_Asymmetry_Results.mat'));

band_names_lower = lower(bandNames);
n_bands = length(bandNames);
n_lobes = length(lobe_labels);
n_subs  = length(unique_subs);

stability_metrics = {'std', 'cv', 'iqr'};
results_cell = {};
pool_idx = struct('gfp_hfp', [], 'ai_lobe_hemi', [], 'ai_channel', [], ...
    'gfp_hfp_stability', [], 'ai_lobe_hemi_stability', [], 'ai_channel_stability', []);

fprintf('\n========== GROUP STATISTICS (Task vs. Rest) ==========\n\n');

% -------------------------------------------------------------------------
% Level tests: GFP/HFP
% -------------------------------------------------------------------------

results_cell{end+1} = run_wilcoxon( ...
    group_gfp_stats.rest.broadband_mean(:), group_gfp_stats.task.broadband_mean(:), ...
    'GFP', 'broadband', '', 'GFP_broadband', 'level', '');

for b = 1:n_bands
    bn = band_names_lower{b};
    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(bn).gfp_mean(:), group_gfp_stats.task.(bn).gfp_mean(:), ...
        'GFP', 'band', bandNames{b}, ['GFP_' bandNames{b}], 'level', '');
end

results_cell{end+1} = run_wilcoxon( ...
    group_gfp_stats.rest.hfp_left_broadband_mean(:), group_gfp_stats.task.hfp_left_broadband_mean(:), ...
    'HFP_Left', 'broadband', '', 'HFP_Left_broadband', 'level', '');

for b = 1:n_bands
    bn = band_names_lower{b};
    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(bn).hfp_left_mean(:), group_gfp_stats.task.(bn).hfp_left_mean(:), ...
        'HFP_Left', 'band', bandNames{b}, ['HFP_Left_' bandNames{b}], 'level', '');
end

results_cell{end+1} = run_wilcoxon( ...
    group_gfp_stats.rest.hfp_right_broadband_mean(:), group_gfp_stats.task.hfp_right_broadband_mean(:), ...
    'HFP_Right', 'broadband', '', 'HFP_Right_broadband', 'level', '');

for b = 1:n_bands
    bn = band_names_lower{b};
    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(bn).hfp_right_mean(:), group_gfp_stats.task.(bn).hfp_right_mean(:), ...
        'HFP_Right', 'band', bandNames{b}, ['HFP_Right_' bandNames{b}], 'level', '');
end

pool_idx.gfp_hfp = 1:length(results_cell);
results_cell = apply_fdr_pool(results_cell, pool_idx.gfp_hfp, 'FDR-BH (GFP+HFP level pool, 18 tests)');

% -------------------------------------------------------------------------
% Level tests: AI lobe/hemi
% -------------------------------------------------------------------------

ai_level_start = length(results_cell) + 1;
for l = 1:n_lobes
    for b = 1:n_bands
        lbl = [lobe_labels{l} '_' bandNames{b}];
        results_cell{end+1} = run_wilcoxon( ...
            squeeze(group_AI_lobe_rest(:, l, b)), squeeze(group_AI_lobe_task(:, l, b)), ...
            'AI_Lobe', lobe_labels{l}, bandNames{b}, lbl, 'level', '');
    end
end

for b = 1:n_bands
    lbl = ['Hemi_' bandNames{b}];
    results_cell{end+1} = run_wilcoxon( ...
        squeeze(group_AI_hemi_rest(:, b)), squeeze(group_AI_hemi_task(:, b)), ...
        'AI_Hemi', 'Hemisphere', bandNames{b}, lbl, 'level', '');
end

pool_idx.ai_lobe_hemi = ai_level_start:length(results_cell);
results_cell = apply_fdr_pool(results_cell, pool_idx.ai_lobe_hemi, 'FDR-BH (AI lobe+hemi level pool, 30 tests)');

% -------------------------------------------------------------------------
% Level tests: AI channel max-T
% -------------------------------------------------------------------------

fprintf('Running max-t paired permutation for AI channel level (%d subjects)...\n', n_subs);
ai_channel_start = length(results_cell) + 1;
[results_cell, pool_idx.ai_channel] = append_channel_max_t(results_cell, ...
    {group_AI_channel_rest}, {group_AI_channel_task}, {''}, ...
    pair_left_labels_ref, bandNames, unique_subs, n_subs, 'level');

% -------------------------------------------------------------------------
% Temporal-stability tests: GFP/HFP
% -------------------------------------------------------------------------

gfp_stab_start = length(results_cell) + 1;
for sm = 1:numel(stability_metrics)
    sm_name = stability_metrics{sm};
    suffix = ['temporal_' sm_name];

    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(['broadband_' suffix])(:), group_gfp_stats.task.(['broadband_' suffix])(:), ...
        'GFP', 'broadband', '', ['GFP_Stability_' sm_name '_broadband'], 'temporal_stability', sm_name);

    for b = 1:n_bands
        bn = band_names_lower{b};
        results_cell{end+1} = run_wilcoxon( ...
            group_gfp_stats.rest.(bn).(['gfp_' suffix])(:), group_gfp_stats.task.(bn).(['gfp_' suffix])(:), ...
            'GFP', 'band', bandNames{b}, ['GFP_Stability_' sm_name '_' bandNames{b}], ...
            'temporal_stability', sm_name);
    end

    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(['hfp_left_broadband_' suffix])(:), group_gfp_stats.task.(['hfp_left_broadband_' suffix])(:), ...
        'HFP_Left', 'broadband', '', ['HFP_Left_Stability_' sm_name '_broadband'], ...
        'temporal_stability', sm_name);

    for b = 1:n_bands
        bn = band_names_lower{b};
        results_cell{end+1} = run_wilcoxon( ...
            group_gfp_stats.rest.(bn).(['hfp_left_' suffix])(:), group_gfp_stats.task.(bn).(['hfp_left_' suffix])(:), ...
            'HFP_Left', 'band', bandNames{b}, ['HFP_Left_Stability_' sm_name '_' bandNames{b}], ...
            'temporal_stability', sm_name);
    end

    results_cell{end+1} = run_wilcoxon( ...
        group_gfp_stats.rest.(['hfp_right_broadband_' suffix])(:), group_gfp_stats.task.(['hfp_right_broadband_' suffix])(:), ...
        'HFP_Right', 'broadband', '', ['HFP_Right_Stability_' sm_name '_broadband'], ...
        'temporal_stability', sm_name);

    for b = 1:n_bands
        bn = band_names_lower{b};
        results_cell{end+1} = run_wilcoxon( ...
            group_gfp_stats.rest.(bn).(['hfp_right_' suffix])(:), group_gfp_stats.task.(bn).(['hfp_right_' suffix])(:), ...
            'HFP_Right', 'band', bandNames{b}, ['HFP_Right_Stability_' sm_name '_' bandNames{b}], ...
            'temporal_stability', sm_name);
    end
end

pool_idx.gfp_hfp_stability = gfp_stab_start:length(results_cell);
results_cell = apply_fdr_pool(results_cell, pool_idx.gfp_hfp_stability, 'FDR-BH (GFP+HFP temporal-stability pool, 54 tests)');

% -------------------------------------------------------------------------
% Temporal-stability tests: AI lobe/hemi
% -------------------------------------------------------------------------

ai_stab_start = length(results_cell) + 1;
for sm = 1:numel(stability_metrics)
    sm_name = stability_metrics{sm};
    for l = 1:n_lobes
        for b = 1:n_bands
            lbl = ['AI_Lobe_Stability_' sm_name '_' lobe_labels{l} '_' bandNames{b}];
            results_cell{end+1} = run_wilcoxon( ...
                squeeze(group_AI_lobe_rest_stability.(sm_name)(:, l, b)), ...
                squeeze(group_AI_lobe_task_stability.(sm_name)(:, l, b)), ...
                'AI_Lobe', lobe_labels{l}, bandNames{b}, lbl, 'temporal_stability', sm_name);
        end
    end

    for b = 1:n_bands
        lbl = ['Hemi_Stability_' sm_name '_' bandNames{b}];
        results_cell{end+1} = run_wilcoxon( ...
            squeeze(group_AI_hemi_rest_stability.(sm_name)(:, b)), ...
            squeeze(group_AI_hemi_task_stability.(sm_name)(:, b)), ...
            'AI_Hemi', 'Hemisphere', bandNames{b}, lbl, 'temporal_stability', sm_name);
    end
end

pool_idx.ai_lobe_hemi_stability = ai_stab_start:length(results_cell);
results_cell = apply_fdr_pool(results_cell, pool_idx.ai_lobe_hemi_stability, 'FDR-BH (AI lobe+hemi temporal-stability pool, 90 tests)');

% -------------------------------------------------------------------------
% Temporal-stability tests: AI channel max-T
% -------------------------------------------------------------------------

fprintf('Running max-t paired permutation for AI channel temporal stability (%d subjects)...\n', n_subs);
[results_cell, pool_idx.ai_channel_stability] = append_channel_max_t(results_cell, ...
    {group_AI_channel_rest_stability.std, group_AI_channel_rest_stability.cv, group_AI_channel_rest_stability.iqr}, ...
    {group_AI_channel_task_stability.std, group_AI_channel_task_stability.cv, group_AI_channel_task_stability.iqr}, ...
    stability_metrics, pair_left_labels_ref, bandNames, unique_subs, n_subs, 'temporal_stability');

% -------------------------------------------------------------------------
% Print and save
% -------------------------------------------------------------------------

fprintf('\n%-20s %-18s %-14s %-9s %-5s %10s %10s %9s %3s %8s %8s %5s\n', ...
    'Type', 'Metric', 'Level', 'Band', 'Stab', 'Mean_Rest', 'Mean_Task', 'Hedges_gz', 'N', 'p_raw', 'p_corr', 'sig');
fprintf('%s\n', repmat('-', 1, 135));

for i = 1:length(results_cell)
    r = results_cell{i};
    sig_mark = '';
    if r.h_fdr && ~isnan(r.p_fdr), sig_mark = '*'; end
    fprintf('%-20s %-18s %-14s %-9s %-5s %10.4f %10.4f %9.3f %3d %8.5f %8.5f %5s\n', ...
        r.analysis_type, r.metric, r.level, r.band, r.stability_metric, ...
        r.mean_rest, r.mean_task, r.effect_size_g, r.n, r.p_raw, r.p_fdr, sig_mark);
end

fprintf('\n%s\n', repmat('-', 1, 135));
fprintf('GFP+HFP level pool:             %d tests (FDR-BH)\n', length(pool_idx.gfp_hfp));
fprintf('AI lobe+hemi level pool:        %d tests (FDR-BH)\n', length(pool_idx.ai_lobe_hemi));
fprintf('AI channel level pool:          %d tests (max-stat permutation)\n', length(pool_idx.ai_channel));
fprintf('GFP+HFP stability pool:         %d tests (FDR-BH)\n', length(pool_idx.gfp_hfp_stability));
fprintf('AI lobe+hemi stability pool:    %d tests (FDR-BH)\n', length(pool_idx.ai_lobe_hemi_stability));
fprintf('AI channel stability pool:      %d tests (max-stat permutation)\n', length(pool_idx.ai_channel_stability));

group_stats = struct();
group_stats.results = results_cell;
group_stats.pool_idx = pool_idx;
group_stats.bandNames = bandNames;
group_stats.lobe_labels = lobe_labels;
group_stats.n_subs = n_subs;
group_stats.unique_subs = unique_subs;
group_stats.stability_metrics = stability_metrics;

save(fullfile(results_dir, 'Group_Statistics.mat'), 'group_stats');
fprintf('\nResults saved to %s\n', fullfile(results_dir, 'Group_Statistics.mat'));
fprintf('\n========== DONE ==========\n');

function results_cell = apply_fdr_pool(results_cell, idx, correction_label)
    pvals = cellfun(@(r) double(r.p_raw), results_cell(idx));
    p_fdr = bh_fdr(pvals);
    for i = 1:length(idx)
        results_cell{idx(i)}.p_fdr = p_fdr(i);
        results_cell{idx(i)}.h_fdr = p_fdr(i) < 0.05;
        results_cell{idx(i)}.correction = correction_label;
    end
end

function [results_cell, pool_indices] = append_channel_max_t(results_cell, rest_sets, task_sets, stability_metrics, ...
        pair_labels, bandNames, unique_subs, n_subs, analysis_type)
    start_idx = length(results_cell) + 1;
    n_sets = numel(rest_sets);
    n_pairs = length(pair_labels);
    n_bands = length(bandNames);

    valid = false(n_subs, 1);
    for s = 1:n_subs
        valid(s) = ~isempty(rest_sets{1}{s}) && ~isempty(task_sets{1}{s});
    end
    n_valid = sum(valid);

    if n_valid < 3
        fprintf('  -> WARNING: fewer than 3 valid subjects for channel max-T. Skipping.\n');
        pool_indices = [];
        return;
    end

    delta_stack = NaN(n_valid, n_pairs, n_bands, n_sets);
    idx_v = 1;
    for s = 1:n_subs
        if valid(s)
            for m = 1:n_sets
                delta_stack(idx_v, :, :, m) = task_sets{m}{s} - rest_sets{m}{s};
            end
            idx_v = idx_v + 1;
        end
    end

    valid_mask = ~isnan(delta_stack);
    n_per_test = squeeze(sum(valid_mask, 1));
    obs_t = NaN(n_pairs, n_bands, n_sets);
    for pair = 1:n_pairs
        for band = 1:n_bands
            for m = 1:n_sets
                nv = n_per_test(pair, band, m);
                if nv >= 3
                    dv = squeeze(delta_stack(:, pair, band, m));
                    dv = dv(~isnan(dv));
                    se = std(dv) / sqrt(nv);
                    if se > 0
                        obs_t(pair, band, m) = mean(dv) / se;
                    end
                end
            end
        end
    end

    n_perms_requested = 50000;
    n_perms_exact = 2^n_valid;
    if n_perms_requested >= n_perms_exact
        n_perms = n_perms_exact;
        flips = dec2bin(0:(n_perms - 1), n_valid) - '0';
        flips(flips == 0) = -1;
        correction_label = sprintf('Exact max-t FWER (2^%d=%d)', n_valid, n_perms);
    else
        rng(1, 'twister');
        n_perms = n_perms_requested;
        flips = randi([0, 1], n_perms, n_valid) * 2 - 1;
        correction_label = sprintf('Monte Carlo max-t FWER (%d random sign flips; seed=1)', n_perms);
    end

    fprintf('  Valid subjects: %d | channel pairs: %d | bands: %d | metrics: %d | %s\n', ...
        n_valid, n_pairs, n_bands, n_sets, correction_label);

    delta_2d = reshape(delta_stack, n_valid, []);
    delta_2d(isnan(delta_2d)) = 0;
    n_per_test_vec = reshape(n_per_test, 1, []);
    sumsq_vec = sum(delta_2d.^2, 1);
    valid_test_vec = n_per_test_vec >= 3;

    max_t_stats = zeros(n_perms, 1);
    progress_step = max(1, floor(n_perms / 10));
    fprintf('  Permutations: 0/%d\n', n_perms);
    perm_timer = tic;
    for p = 1:n_perms
        perm_sum = flips(p, :) * delta_2d;
        perm_mean = perm_sum ./ n_per_test_vec;
        perm_var = (sumsq_vec - (perm_sum.^2 ./ n_per_test_vec)) ./ (n_per_test_vec - 1);
        perm_var = max(perm_var, 0);
        perm_se = sqrt(perm_var ./ n_per_test_vec);
        perm_t = perm_mean ./ perm_se;
        perm_t(~valid_test_vec | perm_se <= 0) = NaN;
        max_t_stats(p) = max(abs(perm_t), [], 'omitnan');
        if mod(p, progress_step) == 0 || p == n_perms
            fprintf('  Permutations: %d/%d\n', p, n_perms);
        end
    end
    fprintf('  Permutations complete.\n');
    fprintf('  Permutation runtime: %.2f seconds.\n', toc(perm_timer));

    for m = 1:n_sets
        sm = stability_metrics{m};
        for pair = 1:n_pairs
            for band = 1:n_bands
                t_abs = abs(obs_t(pair, band, m));
                if isnan(t_abs)
                    p_corr = NaN;
                else
                    p_corr = (sum(max_t_stats >= t_abs) + 1) / (n_perms + 1);
                end

                rest_all = NaN(n_subs, 1);
                task_all = NaN(n_subs, 1);
                for s = 1:n_subs
                    if ~isempty(rest_sets{m}{s}) && ~isempty(task_sets{m}{s})
                        rest_all(s) = rest_sets{m}{s}(pair, band);
                        task_all(s) = task_sets{m}{s}(pair, band);
                    end
                end

                if strcmp(analysis_type, 'level')
                    label = [pair_labels{pair} '_' bandNames{band}];
                else
                    label = [pair_labels{pair} '_Stability_' sm '_' bandNames{band}];
                end

                r = make_result(rest_all, task_all, 'AI_Channel', pair_labels{pair}, bandNames{band}, ...
                    label, 'Permutation (max-t paired)', analysis_type, sm);
                r.p_raw = NaN;
                r.p_fdr = p_corr;
                r.h_fdr = p_corr < 0.05;
                r.correction = correction_label;
                results_cell{end+1} = r; %#ok<AGROW>
            end
        end
    end

    pool_indices = start_idx:length(results_cell);
end

function res = run_wilcoxon(rest_vec, task_vec, metric, level, band, label, analysis_type, stability_metric)
    res = make_result(rest_vec, task_vec, metric, level, band, label, 'Wilcoxon', analysis_type, stability_metric);
    if res.n < 3
        return;
    end
    v = ~isnan(double(rest_vec(:))) & ~isnan(double(task_vec(:)));
    [p, ~] = signrank(double(task_vec(v)), double(rest_vec(v)));
    res.p_raw = double(p);
end

function res = make_result(rest_vec, task_vec, metric, level, band, label, test_name, analysis_type, stability_metric)
    rest_vec = double(rest_vec(:));
    task_vec = double(task_vec(:));
    v = ~isnan(rest_vec) & ~isnan(task_vec);
    diff_scores = task_vec(v) - rest_vec(v);
    n = sum(v);

    if n >= 1
        mean_rest = mean(rest_vec(v));
        mean_task = mean(task_vec(v));
        sd_rest = std(rest_vec(v));
        sd_task = std(task_vec(v));
        sd_delta = std(diff_scores);
        iqr_delta = iqr(diff_scores);
    else
        mean_rest = NaN; mean_task = NaN; sd_rest = NaN; sd_task = NaN; sd_delta = NaN; iqr_delta = NaN;
    end

    if n >= 3
        df = n - 1;
        sd_diff = std(diff_scores);
        if sd_diff > 0 && df > 0
            d_z = mean(diff_scores) / sd_diff;
            hedges_gz = d_z * (1 - (3 / (4 * df - 1)));
        else
            hedges_gz = NaN;
        end
        if std(rest_vec(v)) > 0 && std(task_vec(v)) > 0
            spearman_rest_task = corr(rest_vec(v), task_vec(v), 'Type', 'Spearman', 'Rows', 'complete');
        else
            spearman_rest_task = NaN;
        end
    else
        hedges_gz = NaN;
        spearman_rest_task = NaN;
    end

    res = struct('metric', metric, 'level', level, 'band', band, 'label', label, ...
        'analysis_type', analysis_type, 'stability_metric', stability_metric, ...
        'test', test_name, 'mean_rest', mean_rest, 'mean_task', mean_task, ...
        'sd_rest', sd_rest, 'sd_task', sd_task, 'sd_delta', sd_delta, ...
        'iqr_delta', iqr_delta, 'spearman_rest_task', spearman_rest_task, ...
        'effect_size_g', hedges_gz, 'n', n, 'p_raw', NaN, 'p_fdr', NaN, ...
        'h_fdr', false, 'correction', '');
end

function adj_p = bh_fdr(pvals)
    pvals = pvals(:);
    m = numel(pvals);
    valid = ~isnan(pvals);
    adj_p = NaN(m, 1);
    if ~any(valid), return; end

    m_valid = sum(valid);
    [p_sort, idx_sort] = sort(pvals(valid), 'ascend');
    adj_p_v = zeros(m_valid, 1);
    for i = m_valid:-1:1
        if i == m_valid
            adj_p_v(idx_sort(i)) = p_sort(i);
        else
            adj_p_v(idx_sort(i)) = min(adj_p_v(idx_sort(i+1)), p_sort(i) * (m_valid / i));
        end
    end
    adj_p(valid) = adj_p_v;
end
