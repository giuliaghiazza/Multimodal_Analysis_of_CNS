%% Sliding-window ROI-to-ROI Pearson connectivity matrices
% This is a script, not a function.
%
% Input .txt files must contain numeric ROI time series:
%   rows    = time points
%   columns = ROIs
%
% For every selected file, the script computes one ROI x ROI Pearson
% correlation matrix for each sliding window.

clearvars;
clc;

%% Settings
window_size = 30;    % number of time points in each window
step_size   = 1;     % shift between consecutive windows, in time points
output_dir  = fullfile(pwd, 'cnn_sliding_window_output');

%% Select input txt files
[file_names, input_path] = uigetfile('*.txt', ...
    'Select ROI time-series .txt files', ...
    'MultiSelect', 'on');

if isequal(file_names, 0)
    disp('No files selected.');
    return;
end

if ischar(file_names)
    file_names = {file_names};
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

n_files = numel(file_names);
results = struct( ...
    'file_name', cell(n_files, 1), ...
    'file_path', cell(n_files, 1), ...
    'data', cell(n_files, 1), ...
    'window_size', cell(n_files, 1), ...
    'step_size', cell(n_files, 1), ...
    'window_start', cell(n_files, 1), ...
    'window_end', cell(n_files, 1), ...
    'pearson_matrices', cell(n_files, 1), ...
    'mean_pearson_matrix', cell(n_files, 1), ...
    'figure_file', cell(n_files, 1));

%% Process each file
for f = 1:n_files
    file_path = fullfile(input_path, file_names{f});
    fprintf('Processing file %d/%d: %s\n', f, n_files, file_names{f});

    data = readmatrix(file_path);
    data = data(~all(isnan(data), 2), ~all(isnan(data), 1));

    [n_timepoints, n_rois] = size(data);

    if n_timepoints == 0 || n_rois == 0
        warning('Skipping %s: no numeric data found.', file_names{f});
        continue;
    end

    if n_timepoints < window_size
        warning(['Skipping %s: only %d time points, but window_size is %d.'], ...
            file_names{f}, n_timepoints, window_size);
        continue;
    end

    window_start = 1:step_size:(n_timepoints - window_size + 1);
    window_end = window_start + window_size - 1;
    n_windows = numel(window_start);

    pearson_matrices = nan(n_rois, n_rois, n_windows);

    for w = 1:n_windows
        window_data = data(window_start(w):window_end(w), :);

        % Pearson correlation between ROI columns.
        pearson_matrices(:, :, w) = corr(window_data, ...
            'Type', 'Pearson', ...
            'Rows', 'pairwise');
    end

    mean_pearson_matrix = mean(pearson_matrices, 3, 'omitnan');

    results(f).file_name = file_names{f};
    results(f).file_path = file_path;
    results(f).data = data;
    results(f).window_size = window_size;
    results(f).step_size = step_size;
    results(f).window_start = window_start;
    results(f).window_end = window_end;
    results(f).pearson_matrices = pearson_matrices;
    results(f).mean_pearson_matrix = mean_pearson_matrix;

    [~, base_name] = fileparts(file_names{f});
    output_file = fullfile(output_dir, [base_name '_sliding_window_pearson.mat']);
    save(output_file, 'data', 'window_size', 'step_size', ...
        'window_start', 'window_end', 'pearson_matrices', ...
        'mean_pearson_matrix');

    fig = figure('Name', [base_name ' mean ROI connectivity'], 'Color', 'w');
    imagesc(mean_pearson_matrix);
    axis image;
    colorbar;
    colormap(parula);
    caxis([-1 1]);
    xlabel('ROI');
    ylabel('ROI');
    title(sprintf('%s mean sliding-window Pearson connectivity', base_name), ...
        'Interpreter', 'none');

    figure_file = fullfile(output_dir, [base_name '_mean_pearson_connectivity.png']);
    saveas(fig, figure_file);
    results(f).figure_file = figure_file;

    fprintf('  Matrix data: %d time points x %d ROIs\n', n_timepoints, n_rois);
    fprintf('  Saved %d ROI-to-ROI Pearson matrices to %s\n', ...
        n_windows, output_file);
    fprintf('  Saved connectivity figure to %s\n', figure_file);
end

save(fullfile(output_dir, 'all_sliding_window_pearson_results.mat'), 'results');
disp('Done.');
