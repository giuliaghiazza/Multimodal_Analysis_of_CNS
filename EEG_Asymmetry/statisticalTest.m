% GROUP STATISTICS — Task vs. Rest: GFP, HFP, AI (multiple scales)
% Wilcoxon:  GFP/HFP broad+bands, AI lobe+hemi  (FDR-BH per pool)
% Max-t paired permutation:  AI channel level (FWER built-in)
%
% Pool 1 (GFP+HFP):         18 tests  -> FDR-BH
% Pool 2 (AI lobe+hemi):    30 tests  -> FDR-BH
% Pool 3 (AI channel):    ~135 tests  -> max-stat permutation (no extra FDR)

clear; clc;

base_dir    = fileparts(mfilename('fullpath'));
results_dir = fullfile(base_dir, 'RESULTS');

load(fullfile(results_dir, 'Group_GFP_HFP_Results.mat'));
load(fullfile(results_dir, 'Group_Asymmetry_Results.mat'));

band_names_lower = lower(bandNames);
n_bands   = length(bandNames);
n_lobes   = length(lobe_labels);
n_subs    = length(unique_subs);

results_cell = {};
pool_idx = struct('gfp_hfp', [], 'ai_lobe_hemi', [], 'ai_channel', []);

fprintf('\n========== GROUP STATISTICS (Task vs. Rest) ==========\n\n');

% =====================================================================
% 1. GFP (broadband + 5 bands)  —  6 tests
% =====================================================================

rest_v = group_gfp_stats.rest.broadband_mean(:);
task_v = group_gfp_stats.task.broadband_mean(:);
results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'GFP', 'broadband', '', 'GFP_broadband');

for b = 1:n_bands
    bn = band_names_lower{b};
    rest_v = group_gfp_stats.rest.(bn).gfp_mean(:);
    task_v = group_gfp_stats.task.(bn).gfp_mean(:);
    results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'GFP', 'band', bandNames{b}, ['GFP_' bandNames{b}]);
end

% =====================================================================
% 2. HFP Left (broadband + 5 bands)  —  6 tests
% =====================================================================

rest_v = group_gfp_stats.rest.hfp_left_broadband_mean(:);
task_v = group_gfp_stats.task.hfp_left_broadband_mean(:);
results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'HFP_Left', 'broadband', '', 'HFP_Left_broadband');

for b = 1:n_bands
    bn = band_names_lower{b};
    rest_v = group_gfp_stats.rest.(bn).hfp_left_mean(:);
    task_v = group_gfp_stats.task.(bn).hfp_left_mean(:);
    results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'HFP_Left', 'band', bandNames{b}, ['HFP_Left_' bandNames{b}]);
end

% =====================================================================
% 3. HFP Right (broadband + 5 bands)  —  6 tests
% =====================================================================

rest_v = group_gfp_stats.rest.hfp_right_broadband_mean(:);
task_v = group_gfp_stats.task.hfp_right_broadband_mean(:);
results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'HFP_Right', 'broadband', '', 'HFP_Right_broadband');

for b = 1:n_bands
    bn = band_names_lower{b};
    rest_v = group_gfp_stats.rest.(bn).hfp_right_mean(:);
    task_v = group_gfp_stats.task.(bn).hfp_right_mean(:);
    results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'HFP_Right', 'band', bandNames{b}, ['HFP_Right_' bandNames{b}]);
end

pool_idx.gfp_hfp = 1:length(results_cell);

% =====================================================================
% FDR-BH for GFP/HFP pool (18 tests)
% =====================================================================

p_gh = cellfun(@(r) r.p_raw, results_cell(pool_idx.gfp_hfp));
p_gh_fdr = bh_fdr(p_gh);
for i = 1:length(pool_idx.gfp_hfp)
    results_cell{pool_idx.gfp_hfp(i)}.p_fdr    = p_gh_fdr(i);
    results_cell{pool_idx.gfp_hfp(i)}.h_fdr    = p_gh_fdr(i) < 0.05;
    results_cell{pool_idx.gfp_hfp(i)}.correction = 'FDR-BH (GFP+HFP pool, 18 tests)';
end

% =====================================================================
% 4. AI Lobe (5 lobes x 5 bands)  —  25 tests
% =====================================================================

ai_start1 = length(results_cell) + 1;

for l = 1:n_lobes
    for b = 1:n_bands
        rest_v = squeeze(group_AI_lobe_rest(:, l, b));
        task_v = squeeze(group_AI_lobe_task(:, l, b));
        lbl = [lobe_labels{l} '_' bandNames{b}];
        results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'AI_Lobe', lobe_labels{l}, bandNames{b}, lbl);
    end
end

% =====================================================================
% 5. AI Hemisphere (5 bands)  —  5 tests
% =====================================================================

ai_start2 = length(results_cell) + 1;

for b = 1:n_bands
    rest_v = squeeze(group_AI_hemi_rest(:, b));
    task_v = squeeze(group_AI_hemi_task(:, b));
    lbl = ['Hemi_' bandNames{b}];
    results_cell{end+1} = run_wilcoxon(rest_v, task_v, 'AI_Hemi', 'Hemisphere', bandNames{b}, lbl);
end

pool_idx.ai_lobe_hemi = ai_start1:length(results_cell);

% =====================================================================
% FDR-BH for AI lobe+hemi pool (30 tests)
% =====================================================================

p_lh = cellfun(@(r) r.p_raw, results_cell(pool_idx.ai_lobe_hemi));
p_lh_fdr = bh_fdr(p_lh);
for i = 1:length(pool_idx.ai_lobe_hemi)
    results_cell{pool_idx.ai_lobe_hemi(i)}.p_fdr    = p_lh_fdr(i);
    results_cell{pool_idx.ai_lobe_hemi(i)}.h_fdr    = p_lh_fdr(i) < 0.05;
    results_cell{pool_idx.ai_lobe_hemi(i)}.correction = 'FDR-BH (AI lobe+hemi pool, 30 tests)';
end

% =====================================================================
% 6. AI Channel — Max-T Paired Permutation  (2^10 = 1024 permutations)
% =====================================================================

fprintf('Running max-t paired permutation for AI channel level (%d subjects)...\n', n_subs);

ai_ch_start = length(results_cell) + 1;
n_pairs = length(pair_left_labels_ref);

valid = true(n_subs, 1);
for s = 1:n_subs
    valid(s) = ~isempty(group_AI_channel_rest{s}) && ~isempty(group_AI_channel_task{s});
end
n_valid = sum(valid);

if n_valid >= 3
    delta_stack = NaN(n_valid, n_pairs, n_bands);
    idx_v = 1;
    for s = 1:n_subs
        if valid(s)
            delta_stack(idx_v, :, :) = group_AI_channel_task{s} - group_AI_channel_rest{s};
            idx_v = idx_v + 1;
        end
    end

    valid_mask = ~isnan(delta_stack);
    n_per_test = squeeze(sum(valid_mask, 1));
    se_mat = NaN(n_pairs, n_bands);
    obs_t  = NaN(n_pairs, n_bands);

    for pair = 1:n_pairs
        for band = 1:n_bands
            nv = n_per_test(pair, band);
            if nv >= 3
                dv = squeeze(delta_stack(:, pair, band));
                dv = dv(~isnan(dv));
                se_mat(pair, band) = std(dv) / sqrt(nv);
                obs_t(pair, band)   = mean(dv) / se_mat(pair, band);
            end
        end
    end

    delta_zeroed = delta_stack;
    delta_zeroed(~valid_mask) = 0;

    n_perms = 2^n_valid;
    flips = dec2bin(0:(n_perms - 1), n_valid) - '0';
    flips(flips == 0) = -1;

    max_t_stats = zeros(n_perms, 1);
    fprintf('  Permutation %4d/%d', 0, n_perms);
    for p = 1:n_perms
        sf = reshape(flips(p, :), n_valid, 1, 1);
        perm_mean = squeeze(sum(delta_zeroed .* sf, 1)) ./ n_per_test;
        perm_t    = perm_mean ./ se_mat;
        max_t_stats(p) = max(abs(perm_t(:)), [], 'omitnan');
        if mod(p, 128) == 0
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bPermutation %4d/%d', p, n_perms);
        end
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... done.\n');

    for pair = 1:n_pairs
        for band = 1:n_bands
            t_abs = abs(obs_t(pair, band));
            if isnan(t_abs)
                p_corr = NaN;
            else
                p_corr = sum(max_t_stats >= t_abs) / n_perms;
            end

            rest_all = NaN(n_subs, 1);
            task_all = NaN(n_subs, 1);
            for s = 1:n_subs
                if valid(s)
                    rest_all(s) = group_AI_channel_rest{s}(pair, band);
                    task_all(s) = group_AI_channel_task{s}(pair, band);
                end
            end
            v_ai = ~isnan(rest_all) & ~isnan(task_all);

            r = struct();
            r.test       = 'Permutation (max-t paired)';
            r.correction = sprintf('Max-t FWER (2^%d=%d, built-in)', n_valid, n_perms);
            r.metric     = 'AI_Channel';
            r.level      = pair_left_labels_ref{pair};
            r.band       = bandNames{band};
            r.label      = [pair_left_labels_ref{pair} '_' bandNames{band}];
            r.mean_rest  = mean(rest_all(v_ai));
            r.mean_task  = mean(task_all(v_ai));
            r.n          = sum(v_ai);
            r.p_raw      = p_corr;
            r.p_fdr      = p_corr;
            r.h_fdr      = p_corr < 0.05;

            results_cell{end+1} = r;
        end
    end
else
    fprintf('  -> WARNING: fewer than 3 valid subjects for AI channel permutation. Skipping.\n');
end

pool_idx.ai_channel = ai_ch_start:length(results_cell);

% =====================================================================
% PRINT SUMMARY TABLE
% =====================================================================

fprintf('\n%-22s %-14s %-9s %10s %10s %3s %8s %8s %5s\n', ...
    'Metric', 'Level', 'Band', 'Mean_Rest', 'Mean_Task', 'N', 'p_raw', 'p_corr', 'sig');
fprintf('%s\n', repmat('-', 1, 110));

for i = 1:length(results_cell)
    r = results_cell{i};
    sig_mark = '';
    if r.h_fdr && ~isnan(r.p_fdr), sig_mark = '*'; end
    fprintf('%-22s %-14s %-9s %10.4f %10.4f %3d %8.5f %8.5f %5s\n', ...
        r.metric, r.level, r.band, ...
        r.mean_rest, r.mean_task, r.n, ...
        r.p_raw, r.p_fdr, sig_mark);
end

fprintf('\n%s\n', repmat('-', 1, 110));
fprintf('GFP+HFP pool:      %d tests (FDR-BH)\n', length(pool_idx.gfp_hfp));
fprintf('AI lobe+hemi pool: %d tests (FDR-BH)\n', length(pool_idx.ai_lobe_hemi));
fprintf('AI channel pool:   %d tests (max-stat permutation, no extra FDR)\n', length(pool_idx.ai_channel));

% =====================================================================
% SAVE RESULTS
% =====================================================================

group_stats = struct();
group_stats.results     = results_cell;
group_stats.pool_idx    = pool_idx;
group_stats.bandNames   = bandNames;
group_stats.lobe_labels = lobe_labels;
group_stats.n_subs      = n_subs;
group_stats.unique_subs = unique_subs;

save(fullfile(results_dir, 'Group_Statistics.mat'), 'group_stats');
fprintf('\nResults saved to %s\n', fullfile(results_dir, 'Group_Statistics.mat'));

fprintf('\n========== DONE ==========\n');

% =====================================================================
% LOCAL FUNCTIONS
% =====================================================================

function res = run_wilcoxon(rest_vec, task_vec, metric, level, band, label)
    v = ~isnan(rest_vec) & ~isnan(task_vec) 
    if sum(v) < 3
        res = struct('metric', metric, 'level', level, 'band', band, ...
            'label', label, 'test', 'Wilcoxon', ...
            'mean_rest', NaN, 'mean_task', NaN, 'n', sum(v), ...
            'p_raw', NaN, 'p_fdr', NaN, 'h_fdr', false, 'correction', '');
        return;
    end
    [p, ~] = signrank(task_vec(v), rest_vec(v));
    res = struct('metric', metric, 'level', level, 'band', band, ...
        'label', label, 'test', 'Wilcoxon', ...
        'mean_rest', mean(rest_vec(v)), 'mean_task', mean(task_vec(v)), ...
        'n', sum(v), 'p_raw', p, 'p_fdr', NaN, 'h_fdr', false, 'correction', '');
end

function adj_p = bh_fdr(pvals)
    pvals = pvals(:);
    m = numel(pvals);
    valid = ~isnan(pvals);
    adj_p = NaN(m, 1);
    if ~any(valid), return; end

    [p_sort, idx_sort] = sort(pvals(valid), 'ascend');
    adj_p_v = zeros(sum(valid), 1);
    nv = sum(valid);
    for i = nv:-1:1
        if i == nv
            adj_p_v(idx_sort(i)) = p_sort(i);
        else
            adj_p_v(idx_sort(i)) = min(adj_p_v(idx_sort(i+1)), p_sort(i) * (m / i));
        end
    end
    adj_p(valid) = adj_p_v;
end
