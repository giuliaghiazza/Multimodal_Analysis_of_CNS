% GROUP-LEVEL EEG GFP/HFP FEATURES EXTRACTION

clear; clc;

base_dir = fileparts(mfilename('fullpath'));
if isempty(which('eeglab'))
    addpath(fullfile(base_dir, '..', 'eeglab2026.0.0'));
end
processed_dir = fullfile(base_dir, 'PROCESSED_DATA');
results_dir   = fullfile(base_dir, 'RESULTS');

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

processed_files = dir(fullfile(processed_dir, '*.set'));
subjects_list = {};
for f = 1:length(processed_files)
    parts = strsplit(processed_files(f).name, '_');
    subjects_list{end+1} = parts{1}; %#ok<SAGROW>
end
unique_subs = unique(subjects_list);
n_subs = length(unique_subs);

bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
bands_cfg = [[1;4], [4;8], [8;13], [13;30], [30;45]];
n_bands = length(bandNames);
window_sec = 10;

group_gfp_stats = struct();
cond_keys = {'rest', 'task'};
loaded_counts = struct('rest', 0, 'task', 0);

for ck = 1:numel(cond_keys)
    cond_key = cond_keys{ck};
    group_gfp_stats.(cond_key).broadband_mean = NaN(1, n_subs);
    group_gfp_stats.(cond_key).broadband_std = NaN(1, n_subs);
    group_gfp_stats.(cond_key).broadband_windows = cell(1, n_subs);
    group_gfp_stats.(cond_key).broadband_temporal_std = NaN(1, n_subs);
    group_gfp_stats.(cond_key).broadband_temporal_cv = NaN(1, n_subs);
    group_gfp_stats.(cond_key).broadband_temporal_iqr = NaN(1, n_subs);

    group_gfp_stats.(cond_key).hfp_left_broadband_mean = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_left_broadband_windows = cell(1, n_subs);
    group_gfp_stats.(cond_key).hfp_left_broadband_temporal_std = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_left_broadband_temporal_cv = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_left_broadband_temporal_iqr = NaN(1, n_subs);

    group_gfp_stats.(cond_key).hfp_right_broadband_mean = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_right_broadband_windows = cell(1, n_subs);
    group_gfp_stats.(cond_key).hfp_right_broadband_temporal_std = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_right_broadband_temporal_cv = NaN(1, n_subs);
    group_gfp_stats.(cond_key).hfp_right_broadband_temporal_iqr = NaN(1, n_subs);

    for b = 1:n_bands
        b_name = lower(bandNames{b});
        group_gfp_stats.(cond_key).(b_name).gfp_mean = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).gfp_windows = cell(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).gfp_temporal_std = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).gfp_temporal_cv = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).gfp_temporal_iqr = NaN(1, n_subs);

        group_gfp_stats.(cond_key).(b_name).hfp_left_mean = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_left_windows = cell(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_left_temporal_std = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_left_temporal_cv = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_left_temporal_iqr = NaN(1, n_subs);

        group_gfp_stats.(cond_key).(b_name).hfp_right_mean = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_right_windows = cell(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_right_temporal_std = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_right_temporal_cv = NaN(1, n_subs);
        group_gfp_stats.(cond_key).(b_name).hfp_right_temporal_iqr = NaN(1, n_subs);
    end
end

fprintf('Starting Group GFP/HFP extraction loops across %d subjects...\n', n_subs);

for s = 1:n_subs
    sub_id = unique_subs{s};

    for cond = 1:2
        if cond == 1, cond_str = 'task-eyesclosed';  cond_key = 'rest';
        else,         cond_str = 'task-memory';  cond_key = 'task';
        end

        target_filename = sprintf('%s_ses-session1_%s_processed.set', sub_id, cond_str);
        if ~isfile(fullfile(processed_dir, target_filename))
            fprintf('  -> WARNING: Missing file %s\n', target_filename);
            continue;
        end

        EEG = pop_loadset('filename', target_filename, 'filepath', processed_dir);
        loaded_counts.(cond_key) = loaded_counts.(cond_key) + 1;

        eeg_data = double(EEG.data);
        win_samples = round(window_sec * EEG.srate);

        V_mean = mean(eeg_data, 1);
        gfp_broadband = sqrt(mean((eeg_data - V_mean).^2, 1));

        C = size(eeg_data, 1);
        left_idx = [];
        right_idx = [];
        for c = 1:C
            num_val = str2double(regexp(EEG.chanlocs(c).labels, '\d+', 'match', 'once'));
            if isnan(num_val)
                continue;
            elseif mod(num_val, 2) == 0
                right_idx(end+1) = c; %#ok<SAGROW>
            else
                left_idx(end+1) = c; %#ok<SAGROW>
            end
        end

        V_left_m = mean(eeg_data(left_idx, :), 1);
        hfp_L_broad = sqrt(mean((eeg_data(left_idx, :) - V_left_m).^2, 1));
        V_right_m = mean(eeg_data(right_idx, :), 1);
        hfp_R_broad = sqrt(mean((eeg_data(right_idx, :) - V_right_m).^2, 1));

        group_gfp_stats.(cond_key).broadband_mean(s) = mean(gfp_broadband);
        group_gfp_stats.(cond_key).broadband_std(s) = std(gfp_broadband);
        group_gfp_stats.(cond_key) = assign_temporal_summary(group_gfp_stats.(cond_key), ...
            'broadband', gfp_broadband, win_samples, s);

        group_gfp_stats.(cond_key).hfp_left_broadband_mean(s) = mean(hfp_L_broad);
        group_gfp_stats.(cond_key) = assign_temporal_summary(group_gfp_stats.(cond_key), ...
            'hfp_left_broadband', hfp_L_broad, win_samples, s);

        group_gfp_stats.(cond_key).hfp_right_broadband_mean(s) = mean(hfp_R_broad);
        group_gfp_stats.(cond_key) = assign_temporal_summary(group_gfp_stats.(cond_key), ...
            'hfp_right_broadband', hfp_R_broad, win_samples, s);

        epochData = eeg_data';
        for b = 1:n_bands
            b_name = lower(bandNames{b});
            epochFilt = bandpass(epochData, bands_cfg(:, b)', EEG.srate);
            filt_data = epochFilt';

            V_filt_m = mean(filt_data, 1);
            gfp_b = sqrt(mean((filt_data - V_filt_m).^2, 1));

            V_fl_m = mean(filt_data(left_idx, :), 1);
            hfp_L_b = sqrt(mean((filt_data(left_idx, :) - V_fl_m).^2, 1));
            V_fr_m = mean(filt_data(right_idx, :), 1);
            hfp_R_b = sqrt(mean((filt_data(right_idx, :) - V_fr_m).^2, 1));

            group_gfp_stats.(cond_key).(b_name).gfp_mean(s) = mean(gfp_b);
            group_gfp_stats.(cond_key).(b_name) = assign_temporal_summary( ...
                group_gfp_stats.(cond_key).(b_name), 'gfp', gfp_b, win_samples, s);

            group_gfp_stats.(cond_key).(b_name).hfp_left_mean(s) = mean(hfp_L_b);
            group_gfp_stats.(cond_key).(b_name) = assign_temporal_summary( ...
                group_gfp_stats.(cond_key).(b_name), 'hfp_left', hfp_L_b, win_samples, s);

            group_gfp_stats.(cond_key).(b_name).hfp_right_mean(s) = mean(hfp_R_b);
            group_gfp_stats.(cond_key).(b_name) = assign_temporal_summary( ...
                group_gfp_stats.(cond_key).(b_name), 'hfp_right', hfp_R_b, win_samples, s);
        end
    end
end

fprintf('Loaded GFP/HFP datasets: rest=%d, task=%d\n', loaded_counts.rest, loaded_counts.task);
if loaded_counts.rest == 0 || loaded_counts.task == 0
    error('GFP/HFP extraction did not load both conditions. Check task labels in PROCESSED_DATA.');
end

save(fullfile(results_dir, 'Group_GFP_HFP_Results.mat'), ...
    'group_gfp_stats', 'unique_subs', 'bandNames', 'window_sec');
fprintf('\nGroup GFP/HFP feature aggregation and temporal stability saved to Results directory.\n');

function out = assign_temporal_summary(out, prefix, signal_vec, win_samples, subj_idx)
    [window_values, temporal_std, temporal_cv, temporal_iqr] = temporal_summary(signal_vec, win_samples);
    out.([prefix '_windows']){subj_idx} = window_values;
    out.([prefix '_temporal_std'])(subj_idx) = temporal_std;
    out.([prefix '_temporal_cv'])(subj_idx) = temporal_cv;
    out.([prefix '_temporal_iqr'])(subj_idx) = temporal_iqr;
end

function [window_values, temporal_std, temporal_cv, temporal_iqr] = temporal_summary(signal_vec, win_samples)
    signal_vec = double(signal_vec(:));
    n_windows = floor(numel(signal_vec) / win_samples);
    if n_windows < 1
        window_values = NaN(1, 0);
        temporal_std = NaN;
        temporal_cv = NaN;
        temporal_iqr = NaN;
        return;
    end

    signal_vec = signal_vec(1:(n_windows * win_samples));
    window_mat = reshape(signal_vec, win_samples, n_windows);
    window_values = mean(window_mat, 1, 'omitnan');
    temporal_std = std(window_values, 0, 'omitnan');
    temporal_mean = mean(window_values, 'omitnan');
    if isnan(temporal_mean) || abs(temporal_mean) < eps
        temporal_cv = NaN;
    else
        temporal_cv = temporal_std / abs(temporal_mean);
    end
    temporal_iqr = iqr(window_values);
end
