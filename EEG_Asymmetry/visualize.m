% =========================================================================
% VISUALIZATION — GFP, HFP, AI: Rest vs Task (Group D)
% =========================================================================
% Run AFTER statisticalTest.m has produced RESULTS/Group_Statistics.mat.
% Requires EEGLAB on path only for the optional channel-level AI topoplots.
%
% Current upstream feature files contain one summary value per
% subject/condition/band (and per lobe/channel pair for AI). Therefore this
% script visualizes subject-level paired variation, not epoch distributions.

clear; clc;

base_dir    = fileparts(mfilename('fullpath'));
if isempty(which('eeglab'))
    addpath(fullfile(base_dir, '..', 'eeglab2026.0.0'));
end
results_dir = fullfile(base_dir, 'RESULTS');
fig_dir     = fullfile(results_dir, 'FIGURES');

if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

load(fullfile(results_dir, 'Group_GFP_HFP_Results.mat'));
load(fullfile(results_dir, 'Group_Asymmetry_Results.mat'));
load(fullfile(results_dir, 'Group_Statistics.mat'));

band_names_lower = lower(bandNames);
n_bands = length(bandNames);
n_lobes = length(lobe_labels);
n_subs  = length(unique_subs);

results_cell = group_stats.results;

if all(isnan(group_gfp_stats.task.broadband_mean(:)))
    error(['GFP/HFP task data are all NaN in Group_GFP_HFP_Results.mat. ' ...
        'Rerun GFP_HFP.m after confirming it uses task-mathematic.']);
end

if all(isnan(group_AI_hemi_task(:)))
    error(['AI task data are all NaN in Group_Asymmetry_Results.mat. ' ...
        'Rerun AI.m after confirming it uses task-mathematic.']);
end

has_gfp_hfp_stability = isfield(group_gfp_stats.rest, 'broadband_temporal_std');
has_ai_stability = exist('group_AI_hemi_task_stability', 'var') && ...
    exist('group_AI_channel_task_stability', 'var') && ...
    exist('group_AI_lobe_task_stability', 'var');
has_temporal_stability = has_gfp_hfp_stability && has_ai_stability;

stability_metrics = {'std', 'cv', 'iqr'};
stability_titles = {'STD', 'CV', 'IQR'};

metric_colors = struct();
metric_colors.rest = [0.18 0.44 0.72];
metric_colors.task = [0.82 0.26 0.22];

lobe_profile_colors = [
    0.00 0.45 0.70;  % Frontal
    0.84 0.37 0.00;  % Central
    0.47 0.67 0.19;  % Temporal
    0.49 0.18 0.56;  % Parietal
    0.30 0.75 0.93]; % Occipital

div_map = diverging_map(256);
seq_map = parula(256);

band_labels_display = [{'Broadband'} bandNames];
band_keys_display   = [{'broadband'} bandNames];
n_all_bands = length(band_labels_display);

fprintf('Generating subject-aware figures from Group_Statistics.mat...\n');
if has_temporal_stability
    fprintf('Temporal stability uses 10s non-overlapping windows.\n');
else
    fprintf('Temporal stability fields not found; stability figures will be skipped.\n');
end

% -------------------------------------------------------------------------
% Build corrected-p lookup and AI significance matrices.
% -------------------------------------------------------------------------

stat_map = containers.Map();
sig_ai_lobe = false(n_lobes, n_bands);
sig_ai_hemi = false(1, n_bands);
sig_ai_channel = containers.Map();

for i = 1:length(results_cell)
    r = results_cell{i};
    stat_map(r.label) = r;

    if ~is_level_result(r)
        continue;
    end

    if strcmp(r.metric, 'AI_Lobe')
        li = find(strcmp(lobe_labels, r.level), 1);
        bi = find(strcmp(bandNames, r.band), 1);
        if ~isempty(li) && ~isempty(bi)
            sig_ai_lobe(li, bi) = r.h_fdr && ~isnan(r.p_fdr);
        end
    elseif strcmp(r.metric, 'AI_Hemi')
        bi = find(strcmp(bandNames, r.band), 1);
        if ~isempty(bi)
            sig_ai_hemi(bi) = r.h_fdr && ~isnan(r.p_fdr);
        end
    elseif strcmp(r.metric, 'AI_Channel')
        sig_ai_channel(r.label) = r.h_fdr && ~isnan(r.p_fdr);
    end
end

% -------------------------------------------------------------------------
% Gather matrices: rows = subjects, columns = bands.
% -------------------------------------------------------------------------

gfp_rest = NaN(n_subs, n_all_bands);
gfp_task = NaN(n_subs, n_all_bands);
hfp_left_rest  = NaN(n_subs, n_all_bands);
hfp_left_task  = NaN(n_subs, n_all_bands);
hfp_right_rest = NaN(n_subs, n_all_bands);
hfp_right_task = NaN(n_subs, n_all_bands);

gfp_rest(:, 1) = vec_or_nan(group_gfp_stats.rest.broadband_mean, n_subs);
gfp_task(:, 1) = vec_or_nan(group_gfp_stats.task.broadband_mean, n_subs);
hfp_left_rest(:, 1)  = vec_or_nan(group_gfp_stats.rest.hfp_left_broadband_mean, n_subs);
hfp_left_task(:, 1)  = vec_or_nan(group_gfp_stats.task.hfp_left_broadband_mean, n_subs);
hfp_right_rest(:, 1) = vec_or_nan(group_gfp_stats.rest.hfp_right_broadband_mean, n_subs);
hfp_right_task(:, 1) = vec_or_nan(group_gfp_stats.task.hfp_right_broadband_mean, n_subs);

for b = 1:n_bands
    bn = band_names_lower{b};
    col = b + 1;

    gfp_rest(:, col) = vec_or_nan(group_gfp_stats.rest.(bn).gfp_mean, n_subs);
    gfp_task(:, col) = vec_or_nan(group_gfp_stats.task.(bn).gfp_mean, n_subs);

    hfp_left_rest(:, col)  = vec_or_nan(group_gfp_stats.rest.(bn).hfp_left_mean, n_subs);
    hfp_left_task(:, col)  = vec_or_nan(group_gfp_stats.task.(bn).hfp_left_mean, n_subs);
    hfp_right_rest(:, col) = vec_or_nan(group_gfp_stats.rest.(bn).hfp_right_mean, n_subs);
    hfp_right_task(:, col) = vec_or_nan(group_gfp_stats.task.(bn).hfp_right_mean, n_subs);
end

ai_hemi_rest = group_AI_hemi_rest;
ai_hemi_task = group_AI_hemi_task;
ai_lobe_delta = group_AI_lobe_task - group_AI_lobe_rest;
ai_hemi_delta = ai_hemi_task - ai_hemi_rest;

% -------------------------------------------------------------------------
% FIGURE 1: GFP, paired subject values + subject deltas.
% -------------------------------------------------------------------------

fig = figure('Name', 'GFP Rest vs Task - subject variation', ...
    'Color', 'w', 'Position', [80 80 1400 720]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_paired_bands(gca, gfp_rest, gfp_task, band_labels_display, ...
    stat_map, metric_label_keys('GFP', band_keys_display), ...
    'GFP: paired subject values', 'Mean GFP (\muV)', metric_colors);

nexttile;
plot_subject_delta_heatmap(gca, gfp_task - gfp_rest, unique_subs, ...
    band_labels_display, 'GFP \Delta(Task - Rest) by subject', '\Delta GFP (\muV)', div_map);

save_figure(fig, fig_dir, 'GFP_subject_variation');

% -------------------------------------------------------------------------
% FIGURE 2: HFP left and right, paired subject values + subject deltas.
% -------------------------------------------------------------------------

fig = figure('Name', 'HFP Rest vs Task - subject variation', ...
    'Color', 'w', 'Position', [100 70 1500 850]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_paired_bands(gca, hfp_left_rest, hfp_left_task, band_labels_display, ...
    stat_map, metric_label_keys('HFP_Left', band_keys_display), ...
    'HFP Left: paired subject values', 'Mean HFP (\muV)', metric_colors);

nexttile;
plot_paired_bands(gca, hfp_right_rest, hfp_right_task, band_labels_display, ...
    stat_map, metric_label_keys('HFP_Right', band_keys_display), ...
    'HFP Right: paired subject values', 'Mean HFP (\muV)', metric_colors);

nexttile;
plot_subject_delta_heatmap(gca, hfp_left_task - hfp_left_rest, unique_subs, ...
    band_labels_display, 'HFP Left \Delta(Task - Rest) by subject', '\Delta HFP (\muV)', div_map);

nexttile;
plot_subject_delta_heatmap(gca, hfp_right_task - hfp_right_rest, unique_subs, ...
    band_labels_display, 'HFP Right \Delta(Task - Rest) by subject', '\Delta HFP (\muV)', div_map);

save_figure(fig, fig_dir, 'HFP_subject_variation');

% -------------------------------------------------------------------------
% FIGURE 3: AI lobe summary and subject-by-lobe variation.
% -------------------------------------------------------------------------

fig = figure('Name', 'AI Lobe Rest vs Task', ...
    'Color', 'w', 'Position', [120 70 1450 850]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
ai_lobe_mean_delta = squeeze(mean(ai_lobe_delta, 1, 'omitnan'));
plot_group_heatmap(gca, ai_lobe_mean_delta, lobe_labels, bandNames, ...
    sig_ai_lobe, 'AI Lobe mean \Delta(Task - Rest)', '\Delta AI', div_map);

nexttile;
ai_lobe_abs = squeeze(mean(abs(ai_lobe_delta), 3, 'omitnan'));
imagesc(ai_lobe_abs);
colormap(gca, seq_map);
cb = colorbar;
ylabel(cb, 'Mean |\Delta AI| across bands');
set(gca, 'XTick', 1:n_lobes, 'XTickLabel', lobe_labels, ...
    'YTick', 1:n_subs, 'YTickLabel', unique_subs, 'TickDir', 'out');
xtickangle(30);
title('AI lobe change magnitude by subject');
xlabel('Lobe');
ylabel('Subject');

nexttile([1 2]);
plot_lobe_delta_profiles(gca, ai_lobe_delta, lobe_labels, bandNames, lobe_profile_colors);

save_figure(fig, fig_dir, 'AI_lobe_subject_variation');

% -------------------------------------------------------------------------
% FIGURE 4: AI hemisphere paired values + subject deltas.
% -------------------------------------------------------------------------

fig = figure('Name', 'AI Hemisphere Rest vs Task', ...
    'Color', 'w', 'Position', [140 90 1400 720]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_paired_bands(gca, ai_hemi_rest, ai_hemi_task, bandNames, ...
    stat_map, metric_label_keys('Hemi', bandNames), ...
    'AI Hemisphere: paired subject values', 'AI = log(Right) - log(Left)', metric_colors);

nexttile;
plot_subject_delta_heatmap(gca, ai_hemi_delta, unique_subs, bandNames, ...
    'AI Hemisphere \Delta(Task - Rest) by subject', '\Delta AI', div_map);

save_figure(fig, fig_dir, 'AI_hemisphere_subject_variation');

% -------------------------------------------------------------------------
% FIGURE 5: Distribution summaries with boxplots and histograms.
% -------------------------------------------------------------------------

fig = figure('Name', 'Task-Rest Delta Distributions', ...
    'Color', 'w', 'Position', [90 40 1550 980]);
tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_delta_boxplot(gca, gfp_task - gfp_rest, band_labels_display, ...
    'GFP \Delta(Task - Rest): subject boxplots', '\Delta GFP (\muV)', div_map);

nexttile;
plot_delta_histogram(gca, gfp_task - gfp_rest, band_labels_display, ...
    'GFP \Delta(Task - Rest): histograms', '\Delta GFP (\muV)', div_map);

nexttile;
plot_delta_boxplot(gca, hfp_left_task - hfp_left_rest, band_labels_display, ...
    'HFP Left \Delta(Task - Rest): subject boxplots', '\Delta HFP (\muV)', div_map);

nexttile;
plot_delta_histogram(gca, hfp_left_task - hfp_left_rest, band_labels_display, ...
    'HFP Left \Delta(Task - Rest): histograms', '\Delta HFP (\muV)', div_map);

nexttile;
plot_delta_boxplot(gca, hfp_right_task - hfp_right_rest, band_labels_display, ...
    'HFP Right \Delta(Task - Rest): subject boxplots', '\Delta HFP (\muV)', div_map);

nexttile;
plot_delta_histogram(gca, hfp_right_task - hfp_right_rest, band_labels_display, ...
    'HFP Right \Delta(Task - Rest): histograms', '\Delta HFP (\muV)', div_map);

nexttile;
plot_delta_boxplot(gca, ai_hemi_delta, bandNames, ...
    'AI Hemisphere \Delta(Task - Rest): subject boxplots', '\Delta AI', div_map);

nexttile;
plot_delta_histogram(gca, ai_hemi_delta, bandNames, ...
    'AI Hemisphere \Delta(Task - Rest): histograms', '\Delta AI', div_map);

save_figure(fig, fig_dir, 'Delta_boxplots_histograms');

% -------------------------------------------------------------------------
% FIGURE 6: AI channel topoplots per band, with max-t significant pairs.
% -------------------------------------------------------------------------

try
    [chanlocs, chan_labels_all, chan_map, topo_diag] = load_ai_topoplot_template( ...
        base_dir, fig_dir, pair_left_labels_ref, sig_ai_channel, bandNames);
    if ~isempty(chanlocs)
        topo_vals_by_band = cell(n_bands, 1);
        topo_sig_by_band = cell(n_bands, 1);
        for b = 1:n_bands
            [vals_for_topo, sig_chan_idx] = channel_delta_topodata( ...
                group_AI_channel_rest, group_AI_channel_task, pair_left_labels_ref, ...
                chan_map, chan_labels_all, n_subs, b, bandNames{b}, sig_ai_channel);
            topo_vals_by_band{b} = fill_topoplot_neutral(vals_for_topo);
            topo_sig_by_band{b} = sig_chan_idx;
        end
        topo_clim = common_symmetric_clim(topo_vals_by_band);

        fig = figure('Name', 'AI Channel Topoplots', ...
            'Color', 'w', 'Position', [50 50 1350 720]);
        tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        for b = 1:n_bands
            nexttile;
            vals_for_topo = topo_vals_by_band{b};
            sig_chan_idx = topo_sig_by_band{b};

            if isempty(sig_chan_idx)
                topoplot(vals_for_topo, chanlocs, 'electrodes', 'on', 'numcontour', 0);
            else
                topoplot(vals_for_topo, chanlocs, 'electrodes', 'on', 'numcontour', 0, ...
                    'emarker2', {sig_chan_idx, 'o', 'k', 10, 1.8});
            end
            sig_text = significant_pair_title(sig_ai_channel, bandNames{b});
            title([bandNames{b} ' \DeltaAI' sig_text], 'FontSize', 11);
            colormap(gca, div_map);
            caxis(topo_clim);

            if b == n_bands
                cbar = colorbar;
                ylabel(cbar, sprintf('\\Delta AI, common scale [%.3g, %.3g]', topo_clim(1), topo_clim(2)));
            end
        end
        sgtitle(['AI pair \Delta(Task - Rest), mirrored on homologous electrodes; ' ...
            'black rings = max-t FWER p < 0.05']);
        save_figure(fig, fig_dir, 'AI_channel_topoplots');
        fprintf('  -> AI topoplot template: %s\n', topo_diag.template_name);
    else
        fprintf('  -> No processed .set found for topoplot template; skipping channel maps.\n');
    end
catch ME
    fprintf('  -> Topoplot failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% FIGURE 7: GFP/HFP temporal stability.
% -------------------------------------------------------------------------

if has_temporal_stability
    fig = figure('Name', 'GFP/HFP Temporal Stability', ...
        'Color', 'w', 'Position', [60 40 1650 980]);
    tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    stability_plot_specs = {
        'GFP', 'GFP temporal stability', 'GFP';
        'HFP_Left', 'HFP Left temporal stability', 'HFP Left';
        'HFP_Right', 'HFP Right temporal stability', 'HFP Right'};

    for row = 1:size(stability_plot_specs, 1)
        for sm = 1:numel(stability_metrics)
            nexttile;
            [rest_stab, task_stab] = gfp_hfp_stability_matrix( ...
                group_gfp_stats, stability_plot_specs{row, 1}, stability_metrics{sm}, ...
                band_names_lower, n_subs);
            stat_keys = stability_metric_label_keys(stability_plot_specs{row, 1}, ...
                stability_metrics{sm}, band_keys_display);
            plot_paired_bands(gca, rest_stab, task_stab, band_labels_display, stat_map, stat_keys, ...
                [stability_plot_specs{row, 2} ' - ' stability_titles{sm}], ...
                [stability_plot_specs{row, 3} ' ' stability_titles{sm}], metric_colors);
        end
    end

    save_figure(fig, fig_dir, 'GFP_HFP_temporal_stability');
end

% -------------------------------------------------------------------------
% FIGURE 8: AI lobe temporal stability.
% -------------------------------------------------------------------------

if has_temporal_stability
    fig = figure('Name', 'AI Lobe Temporal Stability', ...
        'Color', 'w', 'Position', [90 70 1500 520]);
    tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for sm = 1:numel(stability_metrics)
        sm_name = stability_metrics{sm};
        nexttile;
        delta_mat = squeeze(mean(group_AI_lobe_task_stability.(sm_name) - ...
            group_AI_lobe_rest_stability.(sm_name), 1, 'omitnan'));
        sig_mat = ai_lobe_stability_sig(stat_map, lobe_labels, bandNames, sm_name);
        plot_group_heatmap(gca, delta_mat, lobe_labels, bandNames, sig_mat, ...
            ['AI Lobe temporal stability \Delta - ' stability_titles{sm}], ...
            ['\Delta ' stability_titles{sm}], div_map);
    end

    save_figure(fig, fig_dir, 'AI_lobe_temporal_stability');
end

% -------------------------------------------------------------------------
% FIGURE 9: AI hemisphere temporal stability.
% -------------------------------------------------------------------------

if has_temporal_stability
    fig = figure('Name', 'AI Hemisphere Temporal Stability', ...
        'Color', 'w', 'Position', [100 40 1500 980]);
    tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for sm = 1:numel(stability_metrics)
        sm_name = stability_metrics{sm};
        rest_stab = group_AI_hemi_rest_stability.(sm_name);
        task_stab = group_AI_hemi_task_stability.(sm_name);
        nexttile;
        plot_paired_bands(gca, rest_stab, task_stab, bandNames, stat_map, ...
            hemi_stability_label_keys(sm_name, bandNames), ...
            ['AI Hemisphere temporal stability - ' stability_titles{sm}], ...
            ['AI ' stability_titles{sm}], metric_colors);
        nexttile;
        plot_subject_delta_heatmap(gca, task_stab - rest_stab, unique_subs, bandNames, ...
            ['AI Hemisphere \Delta stability by subject - ' stability_titles{sm}], ...
            ['\Delta ' stability_titles{sm}], div_map);
    end

    save_figure(fig, fig_dir, 'AI_hemisphere_temporal_stability');
end

% -------------------------------------------------------------------------
% FIGURE 10: AI channel temporal-stability topoplots.
% -------------------------------------------------------------------------

try
    if has_temporal_stability
        [chanlocs, chan_labels_all, chan_map] = load_ai_topoplot_template( ...
            base_dir, fig_dir, pair_left_labels_ref, sig_ai_channel, bandNames);
    else
        chanlocs = [];
    end
    if ~isempty(chanlocs)
        fig = figure('Name', 'AI Channel Temporal Stability Topoplots', ...
            'Color', 'w', 'Position', [40 40 1650 980]);
        tiledlayout(fig, numel(stability_metrics), n_bands, 'TileSpacing', 'compact', 'Padding', 'compact');

        for sm = 1:numel(stability_metrics)
            sm_name = stability_metrics{sm};
            topo_vals_by_band = cell(n_bands, 1);
            topo_sig_by_band = cell(n_bands, 1);
            for b = 1:n_bands
                [vals_for_topo, sig_chan_idx] = channel_stability_delta_topodata( ...
                    group_AI_channel_rest_stability.(sm_name), ...
                    group_AI_channel_task_stability.(sm_name), ...
                    pair_left_labels_ref, chan_map, chan_labels_all, n_subs, ...
                    b, bandNames{b}, sm_name, stat_map);
                topo_vals_by_band{b} = fill_topoplot_neutral(vals_for_topo);
                topo_sig_by_band{b} = sig_chan_idx;
            end
            topo_clim = common_symmetric_clim(topo_vals_by_band);

            for b = 1:n_bands
                nexttile;
                vals_for_topo = topo_vals_by_band{b};
                sig_chan_idx = topo_sig_by_band{b};

                if isempty(sig_chan_idx)
                    topoplot(vals_for_topo, chanlocs, 'electrodes', 'on', 'numcontour', 0);
                else
                    topoplot(vals_for_topo, chanlocs, 'electrodes', 'on', 'numcontour', 0, ...
                        'emarker2', {sig_chan_idx, 'o', 'k', 10, 1.8});
                end
                title([stability_titles{sm} ' ' bandNames{b}], 'FontSize', 10);
                colormap(gca, div_map);
                caxis(topo_clim);
                if b == n_bands
                    cbar = colorbar;
                    ylabel(cbar, sprintf('\\Delta %s, common row scale [%.3g, %.3g]', ...
                        stability_titles{sm}, topo_clim(1), topo_clim(2)));
                end
            end
        end
        sgtitle(['AI pair temporal-stability \Delta(Task - Rest), mirrored on homologous electrodes; ' ...
            'black rings = max-t FWER p < 0.05']);
        save_figure(fig, fig_dir, 'AI_channel_temporal_stability_topoplots');
    end
catch ME
    fprintf('  -> Stability topoplot failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Print and save concise significance summary.
% -------------------------------------------------------------------------

summary_file = fullfile(fig_dir, 'significant_results_summary.txt');
fid = fopen(summary_file, 'w');
fprintf(fid, 'Significant corrected results (p_corr < 0.05)\n');
fprintf(fid, 'Generated by visualize.m\n\n');
fprintf('\nSignificant corrected results:\n');

n_sig = 0;
for i = 1:length(results_cell)
    r = results_cell{i};
    if r.h_fdr && ~isnan(r.p_fdr)
        n_sig = n_sig + 1;
        line = sprintf('%-20s %-12s %-12s %-9s %-5s Rest=% .4f Task=% .4f p_corr=%.5f', ...
            result_field(r, 'analysis_type', 'level'), r.metric, r.level, empty_to_dash(r.band), ...
            empty_to_dash(result_field(r, 'stability_metric', '')), r.mean_rest, r.mean_task, r.p_fdr);
        fprintf('  %s\n', line);
        fprintf(fid, '%s\n', line);
    end
end

if n_sig == 0
    fprintf('  None after correction.\n');
    fprintf(fid, 'None after correction.\n');
end

fclose(fid);

fprintf('\nVisualization complete. Figures saved in:\n  %s\n', fig_dir);
fprintf('Summary saved to:\n  %s\n\n', summary_file);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function y = vec_or_nan(x, n)
    y = NaN(n, 1);
    x = x(:);
    n_copy = min(n, numel(x));
    y(1:n_copy) = x(1:n_copy);
end

function keys = metric_label_keys(prefix, band_keys)
    keys = cell(size(band_keys));
    for k = 1:numel(band_keys)
        keys{k} = [prefix '_' band_keys{k}];
    end
end

function plot_paired_bands(ax, rest_mat, task_mat, band_labels, stat_map, stat_keys, ttl, ylbl, colors)
    axes(ax);
    cla(ax);
    hold(ax, 'on');

    n_bands_local = size(rest_mat, 2);
    x_rest = 1:n_bands_local;
    x_task = x_rest + 0.34;
    jitter = linspace(-0.08, 0.08, size(rest_mat, 1))';

    for b = 1:n_bands_local
        valid = ~isnan(rest_mat(:, b)) & ~isnan(task_mat(:, b));
        for s = find(valid)'
            plot(ax, [x_rest(b) x_task(b)] + jitter(s), ...
                [rest_mat(s, b) task_mat(s, b)], ...
                '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);
        end
        scatter(ax, x_rest(b) + jitter(valid), rest_mat(valid, b), 26, ...
            'MarkerFaceColor', colors.rest, 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.78);
        scatter(ax, x_task(b) + jitter(valid), task_mat(valid, b), 26, ...
            'MarkerFaceColor', colors.task, 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.78);

        mean_rest = mean(rest_mat(valid, b), 'omitnan');
        mean_task = mean(task_mat(valid, b), 'omitnan');
        sem_rest = std(rest_mat(valid, b), 0, 'omitnan') / sqrt(sum(valid));
        sem_task = std(task_mat(valid, b), 0, 'omitnan') / sqrt(sum(valid));

        errorbar(ax, x_rest(b), mean_rest, sem_rest, 'o', ...
            'Color', 'k', 'MarkerFaceColor', colors.rest, 'MarkerSize', 7, 'LineWidth', 1.4);
        errorbar(ax, x_task(b), mean_task, sem_task, 'o', ...
            'Color', 'k', 'MarkerFaceColor', colors.task, 'MarkerSize', 7, 'LineWidth', 1.4);
    end

    set(ax, 'XTick', x_rest + 0.17, 'XTickLabel', band_labels, 'TickDir', 'out');
    xlim(ax, [0.55 n_bands_local + 0.8]);
    ylabel(ax, ylbl);
    title(ax, ttl);
    grid(ax, 'on');
    box(ax, 'off');

    yl = ylim(ax);
    y_span = diff(yl);
    if y_span == 0, y_span = max(abs(yl)); end
    if y_span == 0, y_span = 1; end
    ylim(ax, [yl(1) yl(2) + 0.18*y_span]);
    yl = ylim(ax);

    for b = 1:n_bands_local
        if isKey(stat_map, stat_keys{b})
            r = stat_map(stat_keys{b});
            label = p_label(r.p_fdr, r.h_fdr);
            text(ax, x_rest(b) + 0.17, yl(2) - 0.06*diff(yl), label, ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', [0.12 0.12 0.12]);
        end
    end

    legend(ax, cond_labels_for_legend(colors), {'Rest subjects', 'Task subjects', 'Mean \pm SEM'}, ...
        'Location', 'best', 'Box', 'off');
    hold(ax, 'off');
end

function h = cond_labels_for_legend(colors)
    hold on;
    h(1) = scatter(NaN, NaN, 26, 'MarkerFaceColor', colors.rest, 'MarkerEdgeColor', 'w');
    h(2) = scatter(NaN, NaN, 26, 'MarkerFaceColor', colors.task, 'MarkerEdgeColor', 'w');
    h(3) = errorbar(NaN, NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 1.4);
end

function plot_subject_delta_heatmap(ax, delta_mat, subj_labels, band_labels, ttl, cb_label, cmap)
    axes(ax);
    imagesc(delta_mat);
    colormap(ax, cmap);
    max_abs = max(abs(delta_mat(:)), [], 'omitnan');
    if isnan(max_abs) || max_abs == 0, max_abs = 0.1; end
    caxis(ax, [-max_abs max_abs]);
    cb = colorbar(ax);
    ylabel(cb, cb_label);
    set(ax, 'XTick', 1:numel(band_labels), 'XTickLabel', band_labels, ...
        'YTick', 1:numel(subj_labels), 'YTickLabel', subj_labels, 'TickDir', 'out');
    title(ax, ttl);
    xlabel(ax, 'Band');
    ylabel(ax, 'Subject');
end

function plot_group_heatmap(ax, data_mat, row_labels, col_labels, sig_mat, ttl, cb_label, cmap)
    axes(ax);
    imagesc(data_mat);
    colormap(ax, cmap);
    max_abs = max(abs(data_mat(:)), [], 'omitnan');
    if isnan(max_abs) || max_abs == 0, max_abs = 0.1; end
    caxis(ax, [-max_abs max_abs]);
    cb = colorbar(ax);
    ylabel(cb, cb_label);
    set(ax, 'XTick', 1:numel(col_labels), 'XTickLabel', col_labels, ...
        'YTick', 1:numel(row_labels), 'YTickLabel', row_labels, 'TickDir', 'out');
    title(ax, ttl);
    xlabel(ax, 'Band');
    ylabel(ax, 'Lobe');

    hold(ax, 'on');
    for r = 1:size(data_mat, 1)
        for c = 1:size(data_mat, 2)
            if sig_mat(r, c)
                text(ax, c, r, '*', 'FontSize', 18, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'Color', 'k');
            end
        end
    end
    hold(ax, 'off');
end

function plot_lobe_delta_profiles(ax, ai_lobe_delta, lobe_labels, band_labels, line_colors)
    axes(ax);
    cla(ax);
    hold(ax, 'on');
    n_lobes_local = size(ai_lobe_delta, 2);
    n_bands_local = size(ai_lobe_delta, 3);
    x = 1:n_bands_local;
    mu_mat = NaN(n_lobes_local, n_bands_local);
    sd_mat = NaN(n_lobes_local, n_bands_local);
    line_handles = gobjects(n_lobes_local, 1);

    for l = 1:n_lobes_local
        mat = squeeze(ai_lobe_delta(:, l, :));
        mu_mat(l, :) = mean(mat, 1, 'omitnan');
        sd_mat(l, :) = std(mat, 0, 1, 'omitnan');
    end

    for b = 1:n_bands_local
        sd_for_sort = sd_mat(:, b);
        sd_for_sort(isnan(sd_for_sort)) = Inf;
        [~, order] = sort(sd_for_sort, 'descend');
        for oi = 1:numel(order)
            l = order(oi);
            if isnan(mu_mat(l, b)) || isnan(sd_mat(l, b))
                continue;
            end
            color_idx = mod(l - 1, size(line_colors, 1)) + 1;
            line_color = line_colors(color_idx, :);
            errorbar(ax, x(b), mu_mat(l, b), sd_mat(l, b), 'LineStyle', 'none', ...
                'Color', line_color, 'LineWidth', 1.15, ...
                'CapSize', 9, 'HandleVisibility', 'off');
        end
    end

    for l = 1:n_lobes_local
        color_idx = mod(l - 1, size(line_colors, 1)) + 1;
        line_color = line_colors(color_idx, :);
        line_handles(l) = plot(ax, x, mu_mat(l, :), '-o', ...
            'LineWidth', 1.8, 'MarkerSize', 5, ...
            'Color', line_color, 'MarkerFaceColor', line_color);
    end

    yline(ax, 0, 'k--');
    set(ax, 'XTick', x, 'XTickLabel', band_labels, 'TickDir', 'out');
    xlim(ax, [0.75 n_bands_local + 0.25]);
    title(ax, 'AI lobe mean \Delta(Task - Rest) with SD');
    xlabel(ax, 'Band');
    ylabel(ax, '\Delta AI');
    grid(ax, 'on');
    box(ax, 'off');
    legend(ax, line_handles, lobe_labels, 'Location', 'eastoutside', 'Box', 'off');
    hold(ax, 'off');
end

function plot_delta_boxplot(ax, delta_mat, band_labels, ttl, ylbl, cmap)
    axes(ax);
    cla(ax);
    hold(ax, 'on');

    n_bands_local = size(delta_mat, 2);
    cidx = round(linspace(35, size(cmap, 1) - 35, n_bands_local));
    data_cell = cell(1, n_bands_local);

    for b = 1:n_bands_local
        vals = delta_mat(:, b);
        data_cell{b} = vals(~isnan(vals));
    end

    boxplot(ax, delta_mat, 'Labels', band_labels, 'Symbol', 'k+', 'Whisker', 1.5);
    boxes = findobj(ax, 'Tag', 'Box');
    for b = 1:min(numel(boxes), n_bands_local)
        patch(get(boxes(b), 'XData'), get(boxes(b), 'YData'), ...
            cmap(cidx(n_bands_local - b + 1), :), ...
            'FaceAlpha', 0.22, 'EdgeColor', 'none');
    end

    for b = 1:n_bands_local
        vals = data_cell{b};
        if isempty(vals), continue; end
        jitter = linspace(-0.08, 0.08, numel(vals))';
        scatter(ax, b + jitter, vals, 24, ...
            'MarkerFaceColor', cmap(cidx(b), :), 'MarkerEdgeColor', 'w', ...
            'MarkerFaceAlpha', 0.85);
    end

    yline(ax, 0, 'k--');
    title(ax, ttl);
    ylabel(ax, ylbl);
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'TickDir', 'out');
    xtickangle(ax, 25);
    hold(ax, 'off');
end

function plot_delta_histogram(ax, delta_mat, band_labels, ttl, xlbl, cmap)
    axes(ax);
    cla(ax);
    hold(ax, 'on');

    n_bands_local = size(delta_mat, 2);
    cidx = round(linspace(35, size(cmap, 1) - 35, n_bands_local));
    all_vals = delta_mat(:);
    all_vals = all_vals(~isnan(all_vals));

    if isempty(all_vals)
        title(ax, ttl);
        xlabel(ax, xlbl);
        ylabel(ax, 'Subject count');
        hold(ax, 'off');
        return;
    end

    max_abs = max(abs(all_vals));
    if max_abs == 0, max_abs = 0.1; end
    edges = linspace(-max_abs, max_abs, 12);

    for b = 1:n_bands_local
        vals = delta_mat(:, b);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        histogram(ax, vals, edges, ...
            'DisplayStyle', 'stairs', 'LineWidth', 1.6, ...
            'EdgeColor', cmap(cidx(b), :));
    end

    xline(ax, 0, 'k--');
    title(ax, ttl);
    xlabel(ax, xlbl);
    ylabel(ax, 'Subject count');
    grid(ax, 'on');
    box(ax, 'off');
    legend(ax, band_labels, 'Location', 'bestoutside', 'Box', 'off');
    hold(ax, 'off');
end

function [chanlocs, chan_labels_all, chan_map, diag_info] = load_ai_topoplot_template(base_dir, fig_dir, pair_left_labels, sig_map, band_names)
    processed_dir = fullfile(base_dir, 'PROCESSED_DATA');
    template_files = dir(fullfile(processed_dir, '*.set'));

    chanlocs = [];
    chan_labels_all = {};
    chan_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    diag_info = struct('template_name', '', 'missing_required', {{}}, 'significant_pairs', {{}});

    if isempty(template_files)
        return;
    end

    required_labels = required_ai_pair_labels(pair_left_labels);
    best_score = -1;
    best_eeg = [];
    best_file = '';
    best_labels = {};

    for f = 1:numel(template_files)
        try
            EEG_tmp = pop_loadset('filename', template_files(f).name, 'filepath', processed_dir);
            EEG_tmp = refresh_standard_1005_locations(EEG_tmp);
            labels = {EEG_tmp.chanlocs.labels};
            score = count_present_labels(labels, required_labels);
            if score > best_score
                best_score = score;
                best_eeg = EEG_tmp;
                best_file = template_files(f).name;
                best_labels = labels;
            end
            if score == numel(required_labels)
                break;
            end
        catch ME
            fprintf('  -> Skipping topoplot template %s: %s\n', template_files(f).name, ME.message);
        end
    end

    if isempty(best_eeg)
        return;
    end

    chanlocs = best_eeg.chanlocs;
    chan_labels_all = best_labels;
    chan_map = build_channel_map(chan_labels_all);

    missing_required = labels_missing(chan_labels_all, required_labels);
    significant_pairs = significant_ai_pairs(sig_map, band_names);
    diag_info.template_name = best_file;
    diag_info.missing_required = missing_required;
    diag_info.significant_pairs = significant_pairs;

    if ~isempty(missing_required)
        fprintf('  -> Topoplot warning: missing AI pair channel labels in template: %s\n', ...
            strjoin(missing_required, ', '));
    end
    write_topoplot_diagnostics(fig_dir, diag_info, chan_labels_all);
end

function vals = fill_topoplot_neutral(vals)
    vals(isnan(vals)) = 0;
end

function clim = common_symmetric_clim(vals_by_panel)
    all_vals = [];
    for i = 1:numel(vals_by_panel)
        vals = vals_by_panel{i};
        all_vals = [all_vals; vals(:)]; %#ok<AGROW>
    end
    all_vals = all_vals(~isnan(all_vals));
    max_abs = max(abs(all_vals));
    if isempty(max_abs) || isnan(max_abs) || max_abs == 0
        max_abs = 0.1;
    end
    clim = [-max_abs max_abs];
end

function EEG_tmp = refresh_standard_1005_locations(EEG_tmp)
    if exist('pop_chanedit', 'file') ~= 2
        return;
    end
    try
        EEG_tmp = pop_chanedit(EEG_tmp, 'lookup', 'standard_1005.elc');
    catch
        % Keep saved channel locations if the standard lookup is not available.
    end
end

function chan_map = build_channel_map(chan_labels)
    chan_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    for ch = 1:numel(chan_labels)
        key = normalize_chan_label(chan_labels{ch});
        if ~isempty(key) && ~isKey(chan_map, key)
            chan_map(key) = ch;
        end
    end
end

function labels = required_ai_pair_labels(pair_left_labels)
    labels = {};
    for p = 1:numel(pair_left_labels)
        left_label = pair_left_labels{p};
        right_label = homologous_right_label(left_label);
        labels{end+1} = left_label; %#ok<AGROW>
        if ~isempty(right_label)
            labels{end+1} = right_label; %#ok<AGROW>
        end
    end
    labels = unique(labels, 'stable');
end

function n_present = count_present_labels(chan_labels, required_labels)
    chan_map = build_channel_map(chan_labels);
    n_present = 0;
    for i = 1:numel(required_labels)
        if isKey(chan_map, normalize_chan_label(required_labels{i}))
            n_present = n_present + 1;
        end
    end
end

function missing = labels_missing(chan_labels, required_labels)
    chan_map = build_channel_map(chan_labels);
    missing = {};
    for i = 1:numel(required_labels)
        if ~isKey(chan_map, normalize_chan_label(required_labels{i}))
            missing{end+1} = required_labels{i}; %#ok<AGROW>
        end
    end
end

function sig_pairs = significant_ai_pairs(sig_map, band_names)
    sig_pairs = {};
    keys = sig_map.keys;
    for i = 1:numel(keys)
        key = keys{i};
        if ~sig_map(key)
            continue;
        end
        for b = 1:numel(band_names)
            suffix = ['_' band_names{b}];
            if length(key) > length(suffix) && strcmp(key(end-length(suffix)+1:end), suffix)
                left_label = key(1:end-length(suffix));
                right_label = homologous_right_label(left_label);
                sig_pairs{end+1} = sprintf('%s-%s %s', left_label, right_label, band_names{b}); %#ok<AGROW>
                break;
            end
        end
    end
end

function txt = significant_pair_title(sig_map, band_name)
    labels = {};
    keys = sig_map.keys;
    suffix = ['_' band_name];
    for i = 1:numel(keys)
        key = keys{i};
        if ~sig_map(key)
            continue;
        end
        if length(key) > length(suffix) && strcmp(key(end-length(suffix)+1:end), suffix)
            left_label = key(1:end-length(suffix));
            right_label = homologous_right_label(left_label);
            labels{end+1} = [left_label '-' right_label]; %#ok<AGROW>
        end
    end
    if isempty(labels)
        txt = '';
    else
        txt = [' | sig: ' strjoin(labels, ', ')];
    end
end

function write_topoplot_diagnostics(fig_dir, diag_info, chan_labels)
    fid = fopen(fullfile(fig_dir, 'AI_channel_topoplot_diagnostics.txt'), 'w');
    if fid < 0
        return;
    end
    fprintf(fid, 'AI channel topoplot diagnostics\n');
    fprintf(fid, 'Template .set: %s\n', diag_info.template_name);
    fprintf(fid, 'Coordinates: refreshed with standard_1005.elc when available\n\n');
    fprintf(fid, 'Significant AI pairs from Group_Statistics.mat:\n');
    if isempty(diag_info.significant_pairs)
        fprintf(fid, '  none\n');
    else
        for i = 1:numel(diag_info.significant_pairs)
            fprintf(fid, '  %s\n', diag_info.significant_pairs{i});
        end
    end
    fprintf(fid, '\nMissing required pair labels in template:\n');
    if isempty(diag_info.missing_required)
        fprintf(fid, '  none\n');
    else
        for i = 1:numel(diag_info.missing_required)
            fprintf(fid, '  %s\n', diag_info.missing_required{i});
        end
    end
    fprintf(fid, '\nTemplate labels:\n  %s\n', strjoin(chan_labels, ', '));
    fclose(fid);
end

function [vals_for_topo, sig_chan_idx] = channel_delta_topodata(rest_cells, task_cells, pair_left_labels, chan_map, chan_labels_all, n_subs, band_idx, band_name, sig_map)
    n_pairs = length(pair_left_labels);
    vals_for_topo = NaN(1, length(chan_labels_all));
    sig_chan_idx = [];

    for p = 1:n_pairs
        vals = [];
        for s = 1:n_subs
            if ~isempty(rest_cells{s}) && ~isempty(task_cells{s})
                rv = rest_cells{s}(p, band_idx);
                tv = task_cells{s}(p, band_idx);
                if ~isnan(rv) && ~isnan(tv)
                    vals(end+1) = tv - rv; %#ok<AGROW>
                end
            end
        end
        if numel(vals) < 3, continue; end

        ai_delta = mean(vals, 'omitnan');
        left_label = pair_left_labels{p};
        right_label = homologous_right_label(left_label);
        left_idx = chan_index(chan_map, left_label);
        right_idx = chan_index(chan_map, right_label);

        if ~isempty(left_idx)
            vals_for_topo(left_idx) = ai_delta;
        end
        if ~isempty(right_idx)
            vals_for_topo(right_idx) = ai_delta;
        end

        stat_key = [left_label '_' band_name];
        if isKey(sig_map, stat_key) && sig_map(stat_key)
            if ~isempty(left_idx)
                sig_chan_idx(end+1) = left_idx; %#ok<AGROW>
            end
            if ~isempty(right_idx)
                sig_chan_idx(end+1) = right_idx; %#ok<AGROW>
            end
        end
    end
    sig_chan_idx = unique(sig_chan_idx);
end

function [rest_mat, task_mat] = gfp_hfp_stability_matrix(group_gfp_stats, metric_name, sm_name, band_names_lower, n_subs)
    n_bands_local = numel(band_names_lower);
    rest_mat = NaN(n_subs, n_bands_local + 1);
    task_mat = NaN(n_subs, n_bands_local + 1);

    switch metric_name
        case 'GFP'
            broad_field = ['broadband_temporal_' sm_name];
            band_prefix = 'gfp';
        case 'HFP_Left'
            broad_field = ['hfp_left_broadband_temporal_' sm_name];
            band_prefix = 'hfp_left';
        case 'HFP_Right'
            broad_field = ['hfp_right_broadband_temporal_' sm_name];
            band_prefix = 'hfp_right';
        otherwise
            error('Unknown GFP/HFP metric: %s', metric_name);
    end

    rest_mat(:, 1) = vec_or_nan(group_gfp_stats.rest.(broad_field), n_subs);
    task_mat(:, 1) = vec_or_nan(group_gfp_stats.task.(broad_field), n_subs);

    for b = 1:n_bands_local
        field_name = [band_prefix '_temporal_' sm_name];
        bn = band_names_lower{b};
        rest_mat(:, b + 1) = vec_or_nan(group_gfp_stats.rest.(bn).(field_name), n_subs);
        task_mat(:, b + 1) = vec_or_nan(group_gfp_stats.task.(bn).(field_name), n_subs);
    end
end

function keys = stability_metric_label_keys(prefix, sm_name, band_keys)
    keys = cell(size(band_keys));
    for k = 1:numel(band_keys)
        if strcmpi(band_keys{k}, 'broadband')
            keys{k} = [prefix '_Stability_' sm_name '_broadband'];
        else
            keys{k} = [prefix '_Stability_' sm_name '_' band_keys{k}];
        end
    end
end

function keys = hemi_stability_label_keys(sm_name, band_names)
    keys = cell(size(band_names));
    for k = 1:numel(band_names)
        keys{k} = ['Hemi_Stability_' sm_name '_' band_names{k}];
    end
end

function sig_mat = ai_lobe_stability_sig(stat_map, lobe_labels, band_names, sm_name)
    sig_mat = false(numel(lobe_labels), numel(band_names));
    for l = 1:numel(lobe_labels)
        for b = 1:numel(band_names)
            key = ['AI_Lobe_Stability_' sm_name '_' lobe_labels{l} '_' band_names{b}];
            if isKey(stat_map, key)
                r = stat_map(key);
                sig_mat(l, b) = r.h_fdr && ~isnan(r.p_fdr);
            end
        end
    end
end

function [vals_for_topo, sig_chan_idx] = channel_stability_delta_topodata(rest_cells, task_cells, pair_left_labels, chan_map, chan_labels_all, n_subs, band_idx, band_name, sm_name, stat_map)
    n_pairs = length(pair_left_labels);
    vals_for_topo = NaN(1, length(chan_labels_all));
    sig_chan_idx = [];

    for p = 1:n_pairs
        vals = [];
        for s = 1:n_subs
            if ~isempty(rest_cells{s}) && ~isempty(task_cells{s})
                rv = rest_cells{s}(p, band_idx);
                tv = task_cells{s}(p, band_idx);
                if ~isnan(rv) && ~isnan(tv)
                    vals(end+1) = tv - rv; %#ok<AGROW>
                end
            end
        end
        if numel(vals) < 3, continue; end

        ai_delta = mean(vals, 'omitnan');
        left_label = pair_left_labels{p};
        right_label = homologous_right_label(left_label);
        left_idx = chan_index(chan_map, left_label);
        right_idx = chan_index(chan_map, right_label);

        if ~isempty(left_idx)
            vals_for_topo(left_idx) = ai_delta;
        end
        if ~isempty(right_idx)
            vals_for_topo(right_idx) = ai_delta;
        end

        stat_key = [left_label '_Stability_' sm_name '_' band_name];
        if isKey(stat_map, stat_key)
            r = stat_map(stat_key);
            if r.h_fdr && ~isnan(r.p_fdr)
                if ~isempty(left_idx)
                    sig_chan_idx(end+1) = left_idx; %#ok<AGROW>
                end
                if ~isempty(right_idx)
                    sig_chan_idx(end+1) = right_idx; %#ok<AGROW>
                end
            end
        end
    end
    sig_chan_idx = unique(sig_chan_idx);
end

function idx = chan_index(chan_map, label)
    idx = [];
    key = normalize_chan_label(label);
    if ~isempty(key) && isKey(chan_map, key)
        idx = chan_map(key);
    end
end

function key = normalize_chan_label(label)
    if isempty(label)
        key = '';
    else
        key = upper(strtrim(char(label)));
    end
end

function right_label = homologous_right_label(left_label)
    tok = regexp(left_label, '^([A-Za-z]+)(\d+)$', 'tokens', 'once');
    if isempty(tok)
        right_label = '';
    else
        right_label = [tok{1} num2str(str2double(tok{2}) + 1)];
    end
end

function tf = is_level_result(r)
    tf = ~isfield(r, 'analysis_type') || strcmp(r.analysis_type, 'level');
end

function val = result_field(r, field_name, default_value)
    if isfield(r, field_name)
        val = r.(field_name);
    else
        val = default_value;
    end
end

function label = p_label(p_corr, is_sig)
    if isnan(p_corr)
        label = 'n/a';
    elseif is_sig
        label = sprintf('* p=%.3g', p_corr);
    else
        label = sprintf('p=%.3g', p_corr);
    end
end

function out = empty_to_dash(str)
    if isempty(str)
        out = '-';
    else
        out = str;
    end
end

function cmap = diverging_map(n)
    if nargin < 1, n = 256; end
    n1 = floor(n/2);
    n2 = n - n1;
    blue = [0.16 0.36 0.68];
    white = [0.97 0.97 0.95];
    red = [0.78 0.22 0.18];
    cmap = [interp_colors(blue, white, n1); interp_colors(white, red, n2)];
end

function cmap = interp_colors(c1, c2, n)
    t = linspace(0, 1, n)';
    cmap = c1 .* (1 - t) + c2 .* t;
end

function save_figure(fig, fig_dir, name)
    png_file = fullfile(fig_dir, [name '.png']);
    fig_file = fullfile(fig_dir, [name '.fig']);
    try
        exportgraphics(fig, png_file, 'Resolution', 300);
    catch
        saveas(fig, png_file);
    end
    savefig(fig, fig_file);
end
