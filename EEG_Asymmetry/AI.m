% GROUP-LEVEL EEG ASYMMETRY ANALYSIS

clear; clc;

base_dir = fileparts(mfilename('fullpath'));
if isempty(which('eeglab'))
    addpath(fullfile(base_dir, '..', 'eeglab2026.0.0'));
end
processed_dir = fullfile(base_dir, 'PROCESSED_DATA');
results_dir   = fullfile(base_dir, 'RESULTS');
task_cond     = 'task-memory';
window_sec    = 10;

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

bands = struct();
bands(1).name = 'Delta';     bands(1).range = [1  4];
bands(2).name = 'Theta';     bands(2).range = [4  8];
bands(3).name = 'Alpha';     bands(3).range = [8  13];
bands(4).name = 'Beta';      bands(4).range = [13 30];
bands(5).name = 'Gamma';     bands(5).range = [30 45];
band_labels = {bands.name};
bandNames = band_labels;
n_bands = numel(bands);

lobe_labels = {'Frontal', 'Central', 'Temporal', 'Parietal', 'Occipital'};
lobe_prefixes = {{'Fp', 'AF', 'F', 'FC'}; {'C'}; {'FT', 'T', 'TP'}; {'CP', 'P'}; {'PO', 'O'}};
n_lobes = numel(lobe_labels);

processed_files = dir(fullfile(processed_dir, '*.set'));
if isempty(processed_files), error('No processed .set datasets found.'); end

subjects_list = {};
for f = 1:length(processed_files)
    parts = strsplit(processed_files(f).name, '_');
    subjects_list{end+1} = parts{1}; %#ok<SAGROW>
end
unique_subs = unique(subjects_list);
n_subs = length(unique_subs);

group_AI_lobe_rest = NaN(n_subs, n_lobes, n_bands);
group_AI_lobe_task = NaN(n_subs, n_lobes, n_bands);
group_AI_hemi_rest = NaN(n_subs, n_bands);
group_AI_hemi_task = NaN(n_subs, n_bands);

group_AI_channel_rest = cell(n_subs, 1);
group_AI_channel_task = cell(n_subs, 1);
group_pair_labels     = cell(n_subs, 1);

group_AI_lobe_rest_stability = init_stability_arrays([n_subs, n_lobes, n_bands]);
group_AI_lobe_task_stability = init_stability_arrays([n_subs, n_lobes, n_bands]);
group_AI_hemi_rest_stability = init_stability_arrays([n_subs, n_bands]);
group_AI_hemi_task_stability = init_stability_arrays([n_subs, n_bands]);
group_AI_channel_rest_stability = init_stability_cells(n_subs);
group_AI_channel_task_stability = init_stability_cells(n_subs);

group_AI_lobe_rest_windows = cell(n_subs, 1);
group_AI_lobe_task_windows = cell(n_subs, 1);
group_AI_hemi_rest_windows = cell(n_subs, 1);
group_AI_hemi_task_windows = cell(n_subs, 1);
group_AI_channel_rest_windows = cell(n_subs, 1);
group_AI_channel_task_windows = cell(n_subs, 1);

loaded_counts = struct('rest', 0, 'task', 0);

fprintf('Starting Group Asymmetry Extraction for %d subjects...\n', n_subs);

for s = 1:n_subs
    sub_id = unique_subs{s};
    fprintf('Processing Group Features -> %s\n', sub_id);

    for cond = 1:2
        if cond == 1
            cond_str = 'task-eyesclosed';
            cond_key = 'rest';
        else
            cond_str = task_cond;
            cond_key = 'task';
        end

        target_filename = sprintf('%s_ses-session1_%s_processed.set', sub_id, cond_str);
        if ~isfile(fullfile(processed_dir, target_filename))
            fprintf('  -> WARNING: Missing file %s\n', target_filename);
            continue;
        end

        EEG = pop_loadset('filename', target_filename, 'filepath', processed_dir);
        loaded_counts.(cond_key) = loaded_counts.(cond_key) + 1;

        chan_labels = {EEG.chanlocs.labels};
        [pair_left_idx, pair_right_idx, pair_left_labels] = find_homologous_pairs(chan_labels);
        n_pairs = numel(pair_left_idx);
        fprintf('  -> %s: %d homologous pairs detected\n', cond_str, n_pairs);

        [spectra, freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'plot', 'off', 'verbose', 'off');
        power_linear = double(10 .^ (spectra / 10));
        band_power_lin = band_power_from_spectra(power_linear, freqs, bands);

        [AI_channel, AI_lobe, AI_hemisphere, pair_lobe_idx] = compute_ai_scales( ...
            band_power_lin, pair_left_idx, pair_right_idx, pair_left_labels, lobe_prefixes, n_lobes);

        [AI_channel_windows, AI_lobe_windows, AI_hemi_windows] = compute_ai_windows( ...
            EEG, bands, pair_left_idx, pair_right_idx, pair_left_labels, pair_lobe_idx, n_lobes, window_sec);

        channel_stability = summarize_window_cube(AI_channel_windows);
        lobe_stability = summarize_window_cube(AI_lobe_windows);
        hemi_stability = summarize_window_cube(AI_hemi_windows);

        if cond == 1
            group_AI_lobe_rest(s, :, :) = AI_lobe;
            group_AI_hemi_rest(s, :) = AI_hemisphere;
            group_AI_channel_rest{s} = AI_channel;
            group_AI_lobe_rest_windows{s} = AI_lobe_windows;
            group_AI_hemi_rest_windows{s} = AI_hemi_windows;
            group_AI_channel_rest_windows{s} = AI_channel_windows;
            group_AI_lobe_rest_stability = assign_stability_array(group_AI_lobe_rest_stability, lobe_stability, s);
            group_AI_hemi_rest_stability = assign_stability_array(group_AI_hemi_rest_stability, hemi_stability, s);
            group_AI_channel_rest_stability = assign_stability_cell(group_AI_channel_rest_stability, channel_stability, s);
        else
            group_AI_lobe_task(s, :, :) = AI_lobe;
            group_AI_hemi_task(s, :) = AI_hemisphere;
            group_AI_channel_task{s} = AI_channel;
            group_AI_lobe_task_windows{s} = AI_lobe_windows;
            group_AI_hemi_task_windows{s} = AI_hemi_windows;
            group_AI_channel_task_windows{s} = AI_channel_windows;
            group_AI_lobe_task_stability = assign_stability_array(group_AI_lobe_task_stability, lobe_stability, s);
            group_AI_hemi_task_stability = assign_stability_array(group_AI_hemi_task_stability, hemi_stability, s);
            group_AI_channel_task_stability = assign_stability_cell(group_AI_channel_task_stability, channel_stability, s);
            group_pair_labels{s} = pair_left_labels;
        end
    end
end

fprintf('Loaded AI datasets: rest=%d, task=%d\n', loaded_counts.rest, loaded_counts.task);
if loaded_counts.rest == 0 || loaded_counts.task == 0
    error('AI extraction did not load both conditions. Check task labels in PROCESSED_DATA.');
end

ref_idx = find(~cellfun(@isempty, group_pair_labels), 1);
if isempty(ref_idx)
    error('No homologous left/right channel pairs were found in any %s dataset.', task_cond);
end
pair_left_labels_ref = group_pair_labels{ref_idx};

for s = 1:n_subs
    if ~isempty(group_pair_labels{s}) && ~isequal(group_pair_labels{s}, pair_left_labels_ref)
        error('Homologous channel pair labels differ for %s. Align channel pairs before group statistics.', unique_subs{s});
    end
end

save(fullfile(results_dir, 'Group_Asymmetry_Results.mat'), ...
    'unique_subs', 'lobe_labels', 'band_labels', 'bandNames', 'pair_left_labels_ref', 'window_sec', ...
    'group_AI_lobe_rest', 'group_AI_lobe_task', ...
    'group_AI_hemi_rest', 'group_AI_hemi_task', ...
    'group_AI_channel_rest', 'group_AI_channel_task', ...
    'group_AI_lobe_rest_windows', 'group_AI_lobe_task_windows', ...
    'group_AI_hemi_rest_windows', 'group_AI_hemi_task_windows', ...
    'group_AI_channel_rest_windows', 'group_AI_channel_task_windows', ...
    'group_AI_lobe_rest_stability', 'group_AI_lobe_task_stability', ...
    'group_AI_hemi_rest_stability', 'group_AI_hemi_task_stability', ...
    'group_AI_channel_rest_stability', 'group_AI_channel_task_stability');

fprintf('\nGroup Asymmetry extraction and temporal stability complete. Matrix variables exported to Results.\n');

function [pair_left_idx, pair_right_idx, pair_left_labels] = find_homologous_pairs(chan_labels)
    pair_left_idx = [];
    pair_right_idx = [];
    pair_left_labels = {};
    for ch = 1:numel(chan_labels)
        label = chan_labels{ch};
        tok = regexp(label, '^([A-Za-z]+)(\d+)$', 'tokens', 'once');
        if isempty(tok), continue; end
        prefix = tok{1};
        num = str2double(tok{2});
        if isnan(num) || mod(num, 2) == 0, continue; end
        right_label = sprintf('%s%d', prefix, num + 1);
        idx_r = find(strcmpi(chan_labels, right_label), 1);
        if ~isempty(idx_r)
            pair_left_idx(end+1) = ch; %#ok<AGROW>
            pair_right_idx(end+1) = idx_r; %#ok<AGROW>
            pair_left_labels{end+1} = label; %#ok<AGROW>
        end
    end
end

function band_power_lin = band_power_from_spectra(power_linear, freqs, bands)
    n_chan = size(power_linear, 1);
    n_bands = numel(bands);
    band_power_lin = NaN(n_chan, n_bands);
    for b = 1:n_bands
        band_idx = freqs >= bands(b).range(1) & freqs < bands(b).range(2);
        band_power_lin(:, b) = mean(power_linear(:, band_idx), 2, 'omitnan');
    end
end

function [AI_channel, AI_lobe, AI_hemisphere, pair_lobe_idx] = compute_ai_scales( ...
        band_power_lin, pair_left_idx, pair_right_idx, pair_left_labels, lobe_prefixes, n_lobes)
    n_pairs = numel(pair_left_idx);
    n_bands = size(band_power_lin, 2);
    AI_channel = NaN(n_pairs, n_bands);
    for i = 1:n_pairs
        p_left  = band_power_lin(pair_left_idx(i), :);
        p_right = band_power_lin(pair_right_idx(i), :);
        AI_channel(i, :) = log(p_right + eps) - log(p_left + eps);
    end

    pair_lobe_idx = NaN(n_pairs, 1);
    for i = 1:n_pairs
        tok = regexp(pair_left_labels{i}, '^([A-Za-z]+)\d+$', 'tokens', 'once');
        prefix = tok{1};
        for l = 1:n_lobes
            if any(strcmpi(prefix, lobe_prefixes{l}))
                pair_lobe_idx(i) = l;
                break;
            end
        end
    end

    AI_lobe = NaN(n_lobes, n_bands);
    for l = 1:n_lobes
        mask = pair_lobe_idx == l;
        if any(mask)
            AI_lobe(l, :) = mean(AI_channel(mask, :), 1, 'omitnan');
        end
    end
    AI_hemisphere = mean(AI_channel, 1, 'omitnan');
end

function [AI_channel_windows, AI_lobe_windows, AI_hemi_windows] = compute_ai_windows( ...
        EEG, bands, pair_left_idx, pair_right_idx, pair_left_labels, pair_lobe_idx, n_lobes, window_sec)
    win_samples = round(window_sec * EEG.srate);
    n_windows = floor(EEG.pnts / win_samples);
    n_pairs = numel(pair_left_idx);
    n_bands = numel(bands);
    AI_channel_windows = NaN(n_windows, n_pairs, n_bands);
    AI_lobe_windows = NaN(n_windows, n_lobes, n_bands);
    AI_hemi_windows = NaN(n_windows, n_bands);

    for w = 1:n_windows
        idx = ((w - 1) * win_samples + 1):(w * win_samples);
        data_w = double(EEG.data(:, idx));
        [spectra_w, freqs_w] = spectopo(data_w, numel(idx), EEG.srate, 'plot', 'off', 'verbose', 'off');
        power_linear_w = double(10 .^ (spectra_w / 10));
        band_power_w = band_power_from_spectra(power_linear_w, freqs_w, bands);

        AI_channel_w = NaN(n_pairs, n_bands);
        for i = 1:n_pairs
            p_left = band_power_w(pair_left_idx(i), :);
            p_right = band_power_w(pair_right_idx(i), :);
            AI_channel_w(i, :) = log(p_right + eps) - log(p_left + eps);
        end

        AI_lobe_w = NaN(n_lobes, n_bands);
        for l = 1:n_lobes
            mask = pair_lobe_idx == l;
            if any(mask)
                AI_lobe_w(l, :) = mean(AI_channel_w(mask, :), 1, 'omitnan');
            end
        end
        AI_hemi_w = mean(AI_channel_w, 1, 'omitnan');

        AI_channel_windows(w, :, :) = AI_channel_w;
        AI_lobe_windows(w, :, :) = AI_lobe_w;
        AI_hemi_windows(w, :) = AI_hemi_w;
    end
end

function s = init_stability_arrays(sz)
    s = struct('std', NaN(sz), 'cv', NaN(sz), 'iqr', NaN(sz));
end

function s = init_stability_cells(n_subs)
    s = struct('std', {cell(n_subs, 1)}, 'cv', {cell(n_subs, 1)}, 'iqr', {cell(n_subs, 1)});
end

function dst = assign_stability_array(dst, src, subj_idx)
    if ndims(dst.std) == 3
        dst.std(subj_idx, :, :) = src.std;
        dst.cv(subj_idx, :, :) = src.cv;
        dst.iqr(subj_idx, :, :) = src.iqr;
    else
        dst.std(subj_idx, :) = src.std;
        dst.cv(subj_idx, :) = src.cv;
        dst.iqr(subj_idx, :) = src.iqr;
    end
end

function dst = assign_stability_cell(dst, src, subj_idx)
    dst.std{subj_idx} = src.std;
    dst.cv{subj_idx} = src.cv;
    dst.iqr{subj_idx} = src.iqr;
end

function stability = summarize_window_cube(window_cube)
    if isempty(window_cube)
        stability = struct('std', NaN, 'cv', NaN, 'iqr', NaN);
        return;
    end
    stability.std = squeeze(std(window_cube, 0, 1, 'omitnan'));
    mean_vals = squeeze(mean(window_cube, 1, 'omitnan'));
    stability.cv = stability.std ./ abs(mean_vals);
    stability.cv(abs(mean_vals) < eps) = NaN;
    stability.iqr = squeeze(iqr(window_cube, 1));
end
