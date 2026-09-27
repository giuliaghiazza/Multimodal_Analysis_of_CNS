close all
clear
clc

% ==========================================================
% TD-NIRS DPF analysis with curve-specific temporal windows
% ==========================================================
% 
% This code is used to calculate the DPF for each optode placement and for 
% each wavelength.
%
% For each selected threshold and for each single TPSF curve, 
% the window is recalculated from that curve while the
% threshold value is kept fixed.
%
% Thresholds tested: 1%, 3%, 5%, 10% of the individual curve peak.
%
% 

% ==============================
%% Parameters and Settings
% ==============================

rho_mm = 30;              % mm, source-detector distance
rho = rho_mm * 1e-3;      % m
n = 1.4;                  % refractive index
c0 = 299792458;           % m/s
dt = 6e-12;               % 6 ps per channel

% Wavelengths
lambda = [670 730 780 830];

% Settings to read data
Settings.nChannels  = 4096;
Settings.nCurvesTot = 20;
Settings.dataType   = 'uint16';

% Temporal axis
t = (0:Settings.nChannels-1)' * dt;
t_ps = t * 1e12;

% Pre-signal baseline interval.
% This must be checked in the diagnostic plot.
baseline_idx = 1:200;

% Variable time-window parameters.
% Each threshold is applied to every individual baseline-corrected curve.
threshold_percent = [1 3 5 10];
threshold_fraction_vec = threshold_percent ./ 100;
window_margin_ch = 50;     % channels added before/after threshold crossing

% ==============================
%% IRF reference
% ==============================

Settings.inputFile = 'DPF_MsAB0000.DAT';

curves_ref = TDNIRS_dataLoad(Settings);   % [4096 x 20]
curves_ref = double(curves_ref);

% Keep original reshape structure: [channels x repetitions x wavelengths]
curves_ref_3D = reshape(curves_ref, Settings.nChannels, 5, 4);

% ==============================
%% All TD-NIRS files and sorting
% ==============================

files = dir('DPF_MmAB*.DAT');
nFiles = numel(files);

if nFiles == 0
    error('No DPF_MmAB*.DAT files were found in the current folder.');
end

% Sort files by final numeric code: DPF_MmAB0001.DAT, ..., DPF_MmAB0018.DAT
fileCodes = nan(nFiles,1);
for f = 1:nFiles
    tok = regexp(files(f).name, '(\d{4})\.DAT$', 'tokens', 'once');
    if isempty(tok)
        error('Unable to extract the numeric code from file: %s', files(f).name);
    end
    fileCodes(f) = str2double(tok{1});
end
[~, sortIdx] = sort(fileCodes);
files = files(sortIdx);
fileCodes = fileCodes(sortIdx);

% Segment ID = final file code, according to the montage information
segmentID = fileCodes;

segmentOptodeLabels = [
    "S1-D1"
    "S1-D2"
    "S1-D3"
    "S2-D1"
    "S2-D3"
    "S2-D4"
    "S3-D2"
    "S3-D3"
    "S4-D3"
    "S4-D4"
    "S5-D5"
    "S5-D6"
    "S5-D7"
    "S6-D5"
    "S6-D7"
    "S7-D6"
    "S7-D7"
    "S8-D7"
];

% ==============================
%% Curves for diagnostic plots
% ==============================

plot_curves_irf = cell(1, numel(lambda));
plot_curves_tissue = cell(1, numel(lambda));
window_start_sum = zeros(numel(lambda), numel(threshold_fraction_vec));
window_end_sum = zeros(numel(lambda), numel(threshold_fraction_vec));
window_count = zeros(numel(lambda), numel(threshold_fraction_vec));

for wl = 1:numel(lambda)
    for rep = 1:5
        y_ref_plot = preprocess_tpsf(curves_ref_3D(:, rep, wl), baseline_idx);
        y_ref_plot = normalize_tpsf(y_ref_plot);

        if any(y_ref_plot > 0)
            plot_curves_irf{wl} = [plot_curves_irf{wl}, y_ref_plot];
        end
    end
end

for f = 1:nFiles
    Settings.inputFile = files(f).name;
    curves_tmp = double(TDNIRS_dataLoad(Settings));
    curves3D_tmp = reshape(curves_tmp, Settings.nChannels, 5, 4);

    for wl = 1:numel(lambda)
        for rep = 1:5
            y_plot = preprocess_tpsf(curves3D_tmp(:, rep, wl), baseline_idx);

            if any(y_plot > 0)
                for th = 1:numel(threshold_fraction_vec)
                    [~, i1, i2] = define_curve_window( ...
                        y_plot, threshold_fraction_vec(th), window_margin_ch);

                    if ~isnan(i1) && ~isnan(i2)
                        window_start_sum(wl, th) = window_start_sum(wl, th) + i1;
                        window_end_sum(wl, th) = window_end_sum(wl, th) + i2;
                        window_count(wl, th) = window_count(wl, th) + 1;
                    end
                end

                y_plot = normalize_tpsf(y_plot);
                plot_curves_tissue{wl} = [plot_curves_tissue{wl}, y_plot];
            end
        end
    end
end

mean_tissue_window_ps = nan(numel(lambda), numel(threshold_fraction_vec), 2);
for wl = 1:numel(lambda)
    for th = 1:numel(threshold_fraction_vec)
        if window_count(wl, th) > 0
            mean_i1 = round(window_start_sum(wl, th) / window_count(wl, th));
            mean_i2 = round(window_end_sum(wl, th) / window_count(wl, th));
            mean_tissue_window_ps(wl, th, 1) = t_ps(mean_i1);
            mean_tissue_window_ps(wl, th, 2) = t_ps(mean_i2);
        end
    end
end

fprintf('\nBaseline interval: channels %d-%d, %.1f-%.1f ps\n', ...
    baseline_idx(1), baseline_idx(end), t_ps(baseline_idx(1)), t_ps(baseline_idx(end)));
fprintf('Variable temporal windows: threshold = %s %% of each individual curve peak\n', ...
    strjoin(string(threshold_percent), ', '));
fprintf('Window margin: %d channels\n\n', window_margin_ch);

% Diagnostic plot for baseline and variable window logic
plot_all_tpsf_with_variable_windows( ...
    t_ps, lambda, plot_curves_irf, plot_curves_tissue, ...
    baseline_idx, threshold_fraction_vec, threshold_percent, mean_tissue_window_ps);

% ==============================
%% DPF calculation for all thresholds
% ==============================

nThresh = numel(threshold_fraction_vec);

DPF_all_reps_by_thr = nan(nFiles, 5, 4, nThresh);   % file x repetition x wavelength x threshold
DPF_mean_wl_by_thr  = nan(nFiles, 4, nThresh);      % file x wavelength x threshold
DPF_std_wl_by_thr   = nan(nFiles, 4, nThresh);      % file x wavelength x threshold

DPF_mean_allplacements_by_thr = nan(nThresh, 4);    % threshold x wavelength
DPF_std_allplacements_by_thr  = nan(nThresh, 4);    % threshold x wavelength
DPF_760_by_thr = nan(nThresh, 1);
DPF_850_by_thr = nan(nThresh, 1);

mean_t_ref_all_by_thr = nan(5, 4, nThresh);
mean_t_ref_wl_by_thr  = nan(nThresh, 4);

for th = 1:nThresh
    threshold_fraction = threshold_fraction_vec(th);

    fprintf('=== Threshold %g%% ===\n', threshold_percent(th));

    % Mean time of flight of IRF reference.
    % The temporal window is recalculated separately for each IRF curve.
    for wl = 1:numel(lambda)
        for rep = 1:5
            y_ref = curves_ref_3D(:, rep, wl);
            y_ref = preprocess_tpsf(y_ref, baseline_idx);

            [winMask, i1, i2] = define_curve_window( ...
                y_ref, threshold_fraction, window_margin_ch);

            mean_t_ref_all_by_thr(rep, wl, th) = compute_first_moment(t, y_ref, winMask);

            if isnan(mean_t_ref_all_by_thr(rep, wl, th))
                fprintf('Invalid IRF: threshold %g%%, wavelength %d nm, rep %d.\n', ...
                    threshold_percent(th), lambda(wl), rep);
            else
                fprintf('IRF lambda %d nm rep %d: window ch %d-%d, %.1f-%.1f ps\n', ...
                    lambda(wl), rep, i1, i2, t_ps(i1), t_ps(i2));
            end
        end

        mean_t_ref_wl_by_thr(th, wl) = mean(mean_t_ref_all_by_thr(:, wl, th), 'omitnan');
    end

    % DPF calculation for all files.
    for f = 1:nFiles
        Settings.inputFile = files(f).name;

        curves = TDNIRS_dataLoad(Settings);   % [4096 x 20]
        curves = double(curves);

        % Keep original reshape structure: [channels x repetitions x wavelengths]
        curves3D = reshape(curves, Settings.nChannels, 5, 4);

        for wl = 1:numel(lambda)
            for rep = 1:5
                y = curves3D(:, rep, wl);
                y = preprocess_tpsf(y, baseline_idx);

                winMask = define_curve_window(y, threshold_fraction, window_margin_ch);
                mean_t = compute_first_moment(t, y, winMask);

                if isnan(mean_t) || isnan(mean_t_ref_wl_by_thr(th, wl))
                    DPF_all_reps_by_thr(f, rep, wl, th) = NaN;
                    continue;
                end

                % DPF with IRF reference correction
                DPF_all_reps_by_thr(f, rep, wl, th) = ...
                    c0 * (mean_t - mean_t_ref_wl_by_thr(th, wl)) / (n * rho);
            end

            % Mean and std over the 5 repetitions
            DPF_mean_wl_by_thr(f, wl, th) = mean(DPF_all_reps_by_thr(f, :, wl, th), 'omitnan');
            DPF_std_wl_by_thr(f, wl, th)  = std(DPF_all_reps_by_thr(f, :, wl, th), 0, 'omitnan');
        end

        fprintf('File %s processed. Segment %d.\n', files(f).name, segmentID(f));
    end

    % Summary values and interpolation/extrapolation for this threshold
    DPF_mean_allplacements_by_thr(th, :) = mean(DPF_mean_wl_by_thr(:, :, th), 1, 'omitnan');
    DPF_std_allplacements_by_thr(th, :)  = std(DPF_mean_wl_by_thr(:, :, th), 0, 1, 'omitnan');

    % Local linear interpolation for 760 nm using the closest wavelengths
    idx760 = ismember(lambda, [730 780]);
    DPF_760_by_thr(th) = interp1( ...
        lambda(idx760), DPF_mean_allplacements_by_thr(th, idx760), 760, 'linear');

    % Linear extrapolation for 850 nm using the last two available wavelengths
    idx850 = ismember(lambda, [780 830]);
    DPF_850_by_thr(th) = interp1( ...
        lambda(idx850), DPF_mean_allplacements_by_thr(th, idx850), 850, 'linear', 'extrap');

    fprintf('Estimated DPF at 760 nm: %.6f\n', DPF_760_by_thr(th));
    fprintf('Estimated DPF at 850 nm: %.6f\n\n', DPF_850_by_thr(th));
end

% ==============================
%% Saving and visualizing results
% ==============================

T_results = build_results_table( ...
    files, fileCodes, segmentID, threshold_percent, DPF_mean_wl_by_thr, DPF_std_wl_by_thr);

T_summary = build_summary_table( ...
    lambda, threshold_percent, ...
    DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr, ...
    DPF_760_by_thr, DPF_850_by_thr);

disp(T_results);
disp(T_summary);

% Save CSV files
writetable(T_results, 'DPF_results_finestra_var.csv');
writetable(T_summary, 'DPF_summary_finestra_var.csv');
disp('Results saved in DPF_results_finestra_var.csv and DPF_summary_finestra_var.csv');

%% DPF heatmaps by file/segment and wavelength
rowLabels = "Segment " + string(segmentID) + ": " + segmentOptodeLabels(segmentID);

idxHeatmap1 = find(threshold_percent == 1, 1);
if isempty(idxHeatmap1)
    idxHeatmap1 = 1;
end

plot_combined_heatmap_and_threshold_comparison( ...
    lambda, threshold_percent, DPF_mean_wl_by_thr(:, :, idxHeatmap1), ...
    DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr, ...
    DPF_760_by_thr, DPF_850_by_thr, rowLabels);

for th = 1:nThresh
    if th == 1
        continue;
    end

    DPF_mean_wl = DPF_mean_wl_by_thr(:, :, th);

    figure;
    h = heatmap(string(lambda) + " nm", rowLabels, DPF_mean_wl, ...
        'CellLabelFormat', '%.3f');
    gcf.Color = 'w';
    h.FontColor = 'k';
    h.Title = sprintf('Mean DPF heatmap - threshold %g%%', threshold_percent(th));
    h.XLabel = 'Wavelength';
    h.YLabel = 'Segment / Source-Detector';
    h.Colormap = jet(256);
    h.ColorLimits = [min(DPF_mean_wl, [], 'all', 'omitnan'), ...
                     max(DPF_mean_wl, [], 'all', 'omitnan')];
end

%% Threshold summary plot
plot_threshold_summary( ...
    lambda, threshold_percent, ...
    DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr);

% Optional: save figure
% exportgraphics(gcf, 'DPF_threshold_summary_finestra_var.png', 'Resolution', 300);

% ==========================================================
%% Local helper functions
% ==========================================================

function y_out = preprocess_tpsf(y_in, baseline_idx)
    % Baseline correction from a pre-signal interval.

    y_out = double(y_in);
    bkg = mean(y_out(baseline_idx), 'omitnan');
    y_out = y_out - bkg;
    y_out(y_out < 0) = 0;
end

function [winMask, i1, i2] = define_curve_window(y, threshold_fraction, window_margin_ch)
    % Define a temporal window from the main TPSF lobe only.
    % Starting from the peak avoids isolated noise samples in the tail
    % extending the window at very low thresholds.

    winMask = false(size(y));
    i1 = NaN;
    i2 = NaN;

    [peak, peakIdx] = max(y);
    if peak <= 0 || isnan(peak)
        return;
    end

    threshold = threshold_fraction * peak;

    leftIdx = peakIdx;
    while leftIdx > 1 && y(leftIdx - 1) > threshold
        leftIdx = leftIdx - 1;
    end

    rightIdx = peakIdx;
    while rightIdx < numel(y) && y(rightIdx + 1) > threshold
        rightIdx = rightIdx + 1;
    end

    i1 = max(1, leftIdx - window_margin_ch);
    i2 = min(numel(y), rightIdx + window_margin_ch);
    winMask(i1:i2) = true;
end

function mt = compute_first_moment(t, y, winMask)
    % First temporal moment inside the selected window.
    % Returns NaN if no photons/counts are available in the selected window.

    if ~any(winMask)
        mt = NaN;
        return;
    end

    denom = sum(y(winMask), 'omitnan');
    if denom <= 0 || isnan(denom)
        mt = NaN;
    else
        mt = sum(t(winMask) .* y(winMask), 'omitnan') / denom;
    end
end

function y_norm = normalize_tpsf(y_in)
    peak = max(y_in);
    if peak > 0
        y_norm = y_in ./ peak;
    else
        y_norm = zeros(size(y_in));
    end
end

function T_results = build_results_table( ...
    files, fileCodes, segmentID, threshold_percent, DPF_mean_wl_by_thr, DPF_std_wl_by_thr)

    nFiles = numel(files);
    nThresh = numel(threshold_percent);
    nRows = nFiles * nThresh;

    File = strings(nRows, 1);
    FileCode = nan(nRows, 1);
    Segment = nan(nRows, 1);
    Threshold_percent = nan(nRows, 1);
    DPF_670_mean = nan(nRows, 1);
    DPF_670_std = nan(nRows, 1);
    DPF_730_mean = nan(nRows, 1);
    DPF_730_std = nan(nRows, 1);
    DPF_780_mean = nan(nRows, 1);
    DPF_780_std = nan(nRows, 1);
    DPF_830_mean = nan(nRows, 1);
    DPF_830_std = nan(nRows, 1);

    row = 0;
    for th = 1:nThresh
        for f = 1:nFiles
            row = row + 1;
            File(row) = string(files(f).name);
            FileCode(row) = fileCodes(f);
            Segment(row) = segmentID(f);
            Threshold_percent(row) = threshold_percent(th);
            DPF_670_mean(row) = DPF_mean_wl_by_thr(f, 1, th);
            DPF_670_std(row) = DPF_std_wl_by_thr(f, 1, th);
            DPF_730_mean(row) = DPF_mean_wl_by_thr(f, 2, th);
            DPF_730_std(row) = DPF_std_wl_by_thr(f, 2, th);
            DPF_780_mean(row) = DPF_mean_wl_by_thr(f, 3, th);
            DPF_780_std(row) = DPF_std_wl_by_thr(f, 3, th);
            DPF_830_mean(row) = DPF_mean_wl_by_thr(f, 4, th);
            DPF_830_std(row) = DPF_std_wl_by_thr(f, 4, th);
        end
    end

    T_results = table( ...
        File, FileCode, Segment, Threshold_percent, ...
        DPF_670_mean, DPF_670_std, ...
        DPF_730_mean, DPF_730_std, ...
        DPF_780_mean, DPF_780_std, ...
        DPF_830_mean, DPF_830_std);
end

function T_summary = build_summary_table( ...
    lambda, threshold_percent, DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr, ...
    DPF_760_by_thr, DPF_850_by_thr)

    nThresh = numel(threshold_percent);
    nLambda = numel(lambda);
    nRows = nThresh * nLambda;

    Threshold_percent = nan(nRows, 1);
    Wavelength_nm = nan(nRows, 1);
    DPF_mean_allPlacements = nan(nRows, 1);
    DPF_std_allPlacements = nan(nRows, 1);
    DPF_760_interp_forThreshold = nan(nRows, 1);
    DPF_850_extrap_forThreshold = nan(nRows, 1);

    row = 0;
    for th = 1:nThresh
        for wl = 1:nLambda
            row = row + 1;
            Threshold_percent(row) = threshold_percent(th);
            Wavelength_nm(row) = lambda(wl);
            DPF_mean_allPlacements(row) = DPF_mean_allplacements_by_thr(th, wl);
            DPF_std_allPlacements(row) = DPF_std_allplacements_by_thr(th, wl);
            DPF_760_interp_forThreshold(row) = DPF_760_by_thr(th);
            DPF_850_extrap_forThreshold(row) = DPF_850_by_thr(th);
        end
    end

    T_summary = table( ...
        Threshold_percent, Wavelength_nm, ...
        DPF_mean_allPlacements, DPF_std_allPlacements, ...
        DPF_760_interp_forThreshold, DPF_850_extrap_forThreshold);
end

function plot_all_tpsf_with_variable_windows( ...
    t_ps, lambda, plot_curves_irf, plot_curves_tissue, ...
    baseline_idx, threshold_fraction_vec, threshold_percent, mean_tissue_window_ps)

    figure;
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'All normalized TPSFs with curve-specific threshold levels');

    colors = lines(numel(threshold_fraction_vec));

    for wl = 1:numel(lambda)
        nexttile;
        hold on;

        for th = 1:numel(threshold_fraction_vec)
            x1 = mean_tissue_window_ps(wl, th, 1);
            x2 = mean_tissue_window_ps(wl, th, 2);

            if isnan(x1) || isnan(x2)
                continue;
            end

            patch([x1 x2 x2 x1], [0 0 1.05 1.05], colors(th, :), ...
                'FaceAlpha', 0.08, ...
                'EdgeColor', 'none', ...
                'DisplayName', sprintf('Mean %g%% window', threshold_percent(th)));
        end

        if ~isempty(plot_curves_tissue{wl})
            p_tissue = plot(t_ps, plot_curves_tissue{wl}, ...
                'Color', [0.55 0.70 0.95], 'LineWidth', 0.6);
            p_tissue(1).DisplayName = 'Tissue curves';
            if numel(p_tissue) > 1
                set(p_tissue(2:end), 'HandleVisibility', 'off');
            end

            plot(t_ps, mean(plot_curves_tissue{wl}, 2, 'omitnan'), ...
                'Color', [0.00 0.30 0.85], 'LineWidth', 1.8, 'DisplayName', 'Mean tissue');
        end

        if ~isempty(plot_curves_irf{wl})
            p_irf = plot(t_ps, plot_curves_irf{wl}, ...
                'Color', [0.95 0.60 0.60], 'LineWidth', 0.8);
            p_irf(1).DisplayName = 'IRF curves';
            if numel(p_irf) > 1
                set(p_irf(2:end), 'HandleVisibility', 'off');
            end

            plot(t_ps, mean(plot_curves_irf{wl}, 2, 'omitnan'), ...
                'Color', [0.80 0.10 0.10], 'LineWidth', 1.8, 'DisplayName', 'Mean IRF');
        end

        xline(t_ps(baseline_idx(1)), '--', 'Baseline start');
        xline(t_ps(baseline_idx(end)), '--', 'Baseline end');

        for th = 1:numel(threshold_fraction_vec)
            yline(threshold_fraction_vec(th), ':', ...
                sprintf('%g%% threshold', threshold_percent(th)), ...
                'Color', colors(th, :), 'LineWidth', 1.0);
        end

        xlabel('Time [ps]');
        ylabel('Normalized TPSF [a.u.]');
        title(sprintf('\\lambda = %d nm', lambda(wl)));
        ylim([0 1.05]);
        grid on;

        if wl == 1
            legend('Location', 'best');
        end
    end
end

function plot_threshold_summary( ...
    lambda, threshold_percent, DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr)

    figure('Name', 'DPF threshold summary', 'Color', 'w');
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'DPF mean and standard deviation by threshold');

    colors = lines(numel(threshold_percent));

    for wl = 1:numel(lambda)
        nexttile;
        hold on;

        x = threshold_percent(:);
        y = DPF_mean_allplacements_by_thr(:, wl);
        e = DPF_std_allplacements_by_thr(:, wl);

        errorbar(x, y, e, 'o-', ...
            'Color', [0.00 0.30 0.85], ...
            'MarkerFaceColor', [0.00 0.30 0.85], ...
            'LineWidth', 1.4, ...
            'CapSize', 8, ...
            'DisplayName', 'Mean +/- std across placements');

        for th = 1:numel(threshold_percent)
            plot(x(th), y(th), 'o', ...
                'Color', colors(th, :), ...
                'MarkerFaceColor', colors(th, :), ...
                'MarkerSize', 6, ...
                'HandleVisibility', 'off');
        end

        xlabel('Threshold [% of curve peak]');
        ylabel('DPF');
        title(sprintf('\\lambda = %d nm', lambda(wl)));
        grid on;

        if wl == 1
            legend('Location', 'best');
        end
    end
end

function plot_combined_heatmap_and_threshold_comparison( ...
    lambda, threshold_percent, DPF_mean_wl_1pct, DPF_mean_allplacements_by_thr, DPF_std_allplacements_by_thr, ...
    DPF_760_by_thr, DPF_850_by_thr, rowLabels)

    figure('Name', 'DPF 1% heatmap and threshold comparison', 'Color', 'w');
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'DPF heatmap at 1% threshold and threshold comparison');

    hHeatmap = heatmap(tl, string(lambda) + " nm", rowLabels, DPF_mean_wl_1pct, ...
        'CellLabelFormat', '%.3f');
    hHeatmap.Layout.Tile = 1;
    hHeatmap.FontColor = 'k';
    hHeatmap.Title = 'Mean DPF heatmap - threshold 1%';
    hHeatmap.XLabel = 'Wavelength';
    hHeatmap.YLabel = 'Segment / Source-Detector';
    hHeatmap.Colormap = jet(256);

    colorLimits = [min(DPF_mean_wl_1pct, [], 'all', 'omitnan'), ...
                   max(DPF_mean_wl_1pct, [], 'all', 'omitnan')];
    if all(isfinite(colorLimits))
        if colorLimits(1) == colorLimits(2)
            colorLimits = colorLimits + [-0.5 0.5];
        end
        hHeatmap.ColorLimits = colorLimits;
    end

    axComparison = nexttile(tl, 2);
    hold(axComparison, 'on');
    grid(axComparison, 'on');

    colors = lines(numel(threshold_percent));

    for th = 1:numel(threshold_percent)
        y = DPF_mean_allplacements_by_thr(th, :);
        e = DPF_std_allplacements_by_thr(th, :);

        errorbar(axComparison, lambda, y, e, 'o-', ...
            'LineWidth', 1.4, ...
            'CapSize', 8, ...
            'Color', colors(th, :), ...
            'MarkerFaceColor', colors(th, :), ...
            'DisplayName', sprintf('Threshold %g%%', threshold_percent(th)));

        y730 = y(lambda == 730);
        y780 = y(lambda == 780);
        y830 = y(lambda == 830);

        plot(axComparison, [730 760 780], ...
            [y730, DPF_760_by_thr(th), y780], ...
            '-', 'Color', colors(th, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');

        plot(axComparison, [780 830 850], ...
            [y780, y830, DPF_850_by_thr(th)], ...
            '--', 'Color', colors(th, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');

        plot(axComparison, 760, DPF_760_by_thr(th), 's', ...
            'Color', colors(th, :), ...
            'MarkerFaceColor', 'w', ...
            'MarkerSize', 6, ...
            'HandleVisibility', 'off');

        plot(axComparison, 850, DPF_850_by_thr(th), 'd', ...
            'Color', colors(th, :), ...
            'MarkerFaceColor', 'w', ...
            'MarkerSize', 6, ...
            'HandleVisibility', 'off');
    end

    xlabel(axComparison, 'Wavelength [nm]');
    ylabel(axComparison, 'DPF');
    title(axComparison, 'Mean DPF by wavelength and threshold');
    legend(axComparison, 'Location', 'best');
end
