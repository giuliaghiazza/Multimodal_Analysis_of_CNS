close all
clear
clc

% ==========================================================
% TD-NIRS DPF analysis with a single common temporal window
% ==========================================================
%
% Main changes compared with the previous version:
%   1) Baseline is estimated from an early pre-signal region instead of
%      subtracting min(y). This avoids lifting the whole TPSF by a random
%      negative noise sample.
%   2) The first temporal moment is computed inside one fixed temporal
%      window, common to IRF, all tissue curves, all files, and all
%      wavelengths.
%   3) Reshape structure is kept unchanged:
%      [nChannels x 5 repetitions x 4 wavelengths].
%   4) Source-detector distance is kept fixed at 30 mm.
%
% IMPORTANT: Check that baseline_idx is before the signal/IRF arrival.
%            If the diagnostic figure shows signal in baseline_idx,
%            reduce or shift the baseline interval accordingly.

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
% This must be checked in the diagnostic plot below.
baseline_idx = 1:200;

% Common time-window parameters.
% The representative curve is normalized before averaging so that weak
% channels/wavelengths still contribute to the common window definition.
threshold_fraction = 0.01;   % 1% of representative peak
window_margin_ch   = 50;     % channels added before/after threshold crossing

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
    error('Nessun file DPF_MmAB*.DAT trovato nella cartella corrente.');
end

% Sort files by final numeric code: DPF_MmAB0001.DAT, ..., DPF_MmAB0018.DAT
fileCodes = nan(nFiles,1);
for f = 1:nFiles
    tok = regexp(files(f).name, '(\d{4})\.DAT$', 'tokens', 'once');
    if isempty(tok)
        error('Impossibile estrarre il codice numerico dal file: %s', files(f).name);
    end
    fileCodes(f) = str2double(tok{1});
end
[~, sortIdx] = sort(fileCodes);
files = files(sortIdx);
fileCodes = fileCodes(sortIdx);

% Segment ID = final file code, according to the montage information
segmentID = fileCodes;

% ==============================
%% Define ONE common temporal window
% ==============================
%
% One single absolute channel window is defined from representative IRF and
% tissue curves. The SAME window is then used for every IRF repetition,
% tissue repetition, file, and wavelength.

representative_sum = zeros(Settings.nChannels, 1);
nRepresentative = 0;
plot_curves_irf = cell(1, numel(lambda));
plot_curves_tissue = cell(1, numel(lambda));

% Add normalized mean IRF curves from all wavelengths
for wl = 1:numel(lambda)
    for rep = 1:5
        y_ref_plot = preprocess_tpsf(curves_ref_3D(:, rep, wl), baseline_idx);
        y_ref_plot = normalize_tpsf(y_ref_plot);

        if any(y_ref_plot > 0)
            plot_curves_irf{wl} = [plot_curves_irf{wl}, y_ref_plot];
        end
    end

    y_ref_mean = mean(curves_ref_3D(:,:,wl), 2, 'omitnan');
    y_ref_mean = preprocess_tpsf(y_ref_mean, baseline_idx);

    if max(y_ref_mean) > 0
        representative_sum = representative_sum + y_ref_mean ./ max(y_ref_mean);
        nRepresentative = nRepresentative + 1;
    end
end

% Add normalized mean tissue curves from all files and wavelengths
for f = 1:nFiles
    Settings.inputFile = files(f).name;
    curves_tmp = double(TDNIRS_dataLoad(Settings));

    % Keep original reshape structure
    curves3D_tmp = reshape(curves_tmp, Settings.nChannels, 5, 4);

    for wl = 1:numel(lambda)
        for rep = 1:5
            y_plot = preprocess_tpsf(curves3D_tmp(:, rep, wl), baseline_idx);
            y_plot = normalize_tpsf(y_plot);

            if any(y_plot > 0)
                plot_curves_tissue{wl} = [plot_curves_tissue{wl}, y_plot];
            end
        end

        y_mean = mean(curves3D_tmp(:,:,wl), 2, 'omitnan');
        y_mean = preprocess_tpsf(y_mean, baseline_idx);

        if max(y_mean) > 0
            representative_sum = representative_sum + y_mean ./ max(y_mean);
            nRepresentative = nRepresentative + 1;
        end
    end
end

if nRepresentative == 0
    error('Impossibile definire la finestra temporale: nessuna curva valida.');
end

representative_curve = representative_sum ./ nRepresentative;

idx = find(representative_curve > threshold_fraction * max(representative_curve));
if isempty(idx)
    error('Nessuna finestra temporale valida trovata. Ridurre threshold_fraction.');
end

i1 = max(1, idx(1) - window_margin_ch);
i2 = min(Settings.nChannels, idx(end) + window_margin_ch);
winMask = false(Settings.nChannels, 1);
winMask(i1:i2) = true;

fprintf('\nCommon temporal window: channels %d-%d, %.1f-%.1f ps\n', ...
    i1, i2, t_ps(i1), t_ps(i2));
fprintf('Baseline interval: channels %d-%d, %.1f-%.1f ps\n\n', ...
    baseline_idx(1), baseline_idx(end), t_ps(baseline_idx(1)), t_ps(baseline_idx(end)));

% Diagnostic plot for baseline and common window
figure;
plot(t_ps, representative_curve, 'LineWidth', 1.2);
hold on;
xline(t_ps(baseline_idx(1)), '--', 'Baseline start');
xline(t_ps(baseline_idx(end)), '--', 'Baseline end');
xline(t_ps(i1), '-', 'Window start');
xline(t_ps(i2), '-', 'Window end');
xlabel('Time [ps]');
ylabel('Normalized representative TPSF [a.u.]');
title('Common time window and baseline interval');
grid on;

plot_all_tpsf_with_window(t_ps, lambda, plot_curves_irf, plot_curves_tissue, i1, i2);

% ==============================
%% Mean time of flight of IRF reference
% ==============================

mean_t_ref_all = nan(5,4);
mean_t_ref_wl  = nan(1,4);

for wl = 1:numel(lambda)
    for rep = 1:5
        y_ref = curves_ref_3D(:, rep, wl);
        y_ref = preprocess_tpsf(y_ref, baseline_idx);

        mean_t_ref_all(rep, wl) = compute_first_moment(t, y_ref, winMask);
    end

    % Mean of the 5 IRF repetitions for each wavelength
    mean_t_ref_wl(wl) = mean(mean_t_ref_all(:, wl), 'omitnan');
end

% ==============================
%% DPF calculation for all files
% ==============================

DPF_all_reps  = nan(nFiles, 5, 4);   % file x repetition x wavelength
DPF_mean_wl   = nan(nFiles, 4);      % file x wavelength
DPF_std_wl    = nan(nFiles, 4);      % file x wavelength

for f = 1:nFiles
    Settings.inputFile = files(f).name;

    % Reading the curves for the current file
    curves = TDNIRS_dataLoad(Settings);   % [4096 x 20]
    curves = double(curves);

    % Keep original reshape structure: [channels x repetitions x wavelengths]
    curves3D = reshape(curves, Settings.nChannels, 5, 4);

    for wl = 1:numel(lambda)
        for rep = 1:5
            y = curves3D(:, rep, wl);
            y = preprocess_tpsf(y, baseline_idx);

            mean_t = compute_first_moment(t, y, winMask);

            if isnan(mean_t) || isnan(mean_t_ref_wl(wl))
                DPF_all_reps(f, rep, wl) = NaN;
                continue;
            end

            % DPF with IRF reference correction
            DPF_all_reps(f, rep, wl) = c0 * (mean_t - mean_t_ref_wl(wl)) / (n * rho);
        end

        % Mean and std over the 5 repetitions
        DPF_mean_wl(f, wl) = mean(DPF_all_reps(f, :, wl), 'omitnan');
        DPF_std_wl(f, wl)  = std(DPF_all_reps(f, :, wl),  0, 'omitnan');
    end

    fprintf('File %s elaborato. Segmento %d.\n', files(f).name, segmentID(f));
end

% ==============================
%% Summary values and interpolation/extrapolation
% ==============================

DPF_mean_allplacements = mean(DPF_mean_wl, 1, 'omitnan');
DPF_std_allplacements  = std(DPF_mean_wl, 0, 1, 'omitnan');

T_allplacements = table( ...
    lambda(:), ...
    DPF_mean_allplacements(:), ...
    DPF_std_allplacements(:), ...
    'VariableNames', {'Wavelength_nm', 'DPF_mean_allPlacements', 'DPF_std_allPlacements'});

% Local linear interpolation for 760 nm using the closest wavelengths
idx760 = ismember(lambda, [730 780]);
DPF_760 = interp1(lambda(idx760), DPF_mean_allplacements(idx760), 760, 'linear');

% Linear extrapolation for 850 nm using the last two available wavelengths
idx850 = ismember(lambda, [780 830]);
DPF_850 = interp1(lambda(idx850), DPF_mean_allplacements(idx850), 850, 'linear', 'extrap');

% ==============================
%% Saving and visualizing results
% ==============================

T = table( ...
    string({files.name})', ...
    fileCodes, ...
    segmentID, ...
    DPF_mean_wl(:,1), DPF_std_wl(:,1), ...
    DPF_mean_wl(:,2), DPF_std_wl(:,2), ...
    DPF_mean_wl(:,3), DPF_std_wl(:,3), ...
    DPF_mean_wl(:,4), DPF_std_wl(:,4), ...
    'VariableNames', { ...
    'File', 'FileCode', 'Segment', ...
    'DPF_670_mean','DPF_670_std', ...
    'DPF_730_mean','DPF_730_std', ...
    'DPF_780_mean','DPF_780_std', ...
    'DPF_830_mean','DPF_830_std'});

disp(T);
disp(T_allplacements);

fprintf('\nDPF stimata a 760 nm (interpolazione lineare locale 730-780): %.6f\n', DPF_760);
fprintf('DPF stimata a 850 nm (estrapolazione lineare 780-830): %.6f\n', DPF_850);

% Save CSV files
writetable(T, 'DPF_results_common_window.csv');
writetable(T_allplacements, 'DPF_summary_common_window.csv');
disp('Results saved in DPF_results_common_window.csv and DPF_summary_common_window.csv');

%% Heatmap DPF per file/segmento e lunghezza d'onda
rowLabels = "Seg " + string(segmentID) + " | " + string({files.name})';

figure;
h = heatmap(string(lambda) + " nm", rowLabels, DPF_mean_wl, ...
    'CellLabelFormat', '%.3f');
gcf.Color = 'w';
h.FontColor = 'k';
h.Title = 'Heatmap DPF medio per segmento e lunghezza d''onda';
h.XLabel = 'Wavelength';
h.YLabel = 'Segment / File';
h.Colormap = jet(256);
h.ColorLimits = [min(DPF_mean_wl, [], 'all', 'omitnan'), ...
                 max(DPF_mean_wl, [], 'all', 'omitnan')];

% Optional: save figure
% exportgraphics(gcf, 'DPF_heatmap_common_window.png', 'Resolution', 300);

% ==========================================================
%% Local helper functions
% ==========================================================

function y_out = preprocess_tpsf(y_in, baseline_idx)
    % Baseline correction from a pre-signal interval.
    % This is preferable to subtracting min(y), which can be dominated by
    % a single noise sample and can bias the temporal moment.

    y_out = double(y_in);
    bkg = mean(y_out(baseline_idx), 'omitnan');
    y_out = y_out - bkg;
    y_out(y_out < 0) = 0;
end

function mt = compute_first_moment(t, y, winMask)
    % First temporal moment inside a fixed window.
    % Returns NaN if no photons/counts are available in the selected window.

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

function plot_all_tpsf_with_window(t_ps, lambda, plot_curves_irf, plot_curves_tissue, i1, i2)
    figure;
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'All normalized TPSFs with common time window');

    x_window = [t_ps(i1) t_ps(i2) t_ps(i2) t_ps(i1)];
    y_window = [0 0 1.05 1.05];

    for wl = 1:numel(lambda)
        nexttile;
        hold on;

        patch(x_window, y_window, [1.00 0.93 0.75], ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Common window');

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

        xline(t_ps(i1), '--', 'Window start', 'LabelVerticalAlignment', 'middle');
        xline(t_ps(i2), '--', 'Window end', 'LabelVerticalAlignment', 'middle');
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
