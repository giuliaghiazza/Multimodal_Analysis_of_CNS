close all
clear 
clc

% ==============================
%    Parameters and Settings 
% ==============================

rho_mm = 30;              % mm, source-detector distance
rho = rho_mm * 1e-3;      % m
n = 1.4;
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

% ==============================
%%      IRF reference 
% ==============================

Settings.inputFile = 'DPF_MsAB0000.DAT';

curves_ref = TDNIRS_dataLoad(Settings);   % [4096 x 20] - (5 rip x 4 lambda, checked using TDNIRS_dataLoad)
curves_ref = double(curves_ref);
%size(curves_ref)
%plot(curves_ref(:,1:5))

% Subdividing the repetitions by lambda [4096 x 5 x 4]
curves_ref_3D = reshape(curves_ref, Settings.nChannels, 5, 4);

% Mean time of flight for each repetition and wavelength
mean_t_ref_all = zeros(5,4);
mean_t_ref_wl  = zeros(1,4);

for wl = 1:4
    for rep = 1:5
        y_ref = curves_ref_3D(:, rep, wl);

        % Baseline correction
        y_ref = y_ref - min(y_ref);
        y_ref(y_ref < 0) = 0;

        % Mean time of flight of the reference
        mean_t_ref_all(rep, wl) = sum(t .* y_ref) / sum(y_ref);
    end

    % mean of the 5 repetitions of the reference for each wavelength
    mean_t_ref_wl(wl) = mean(mean_t_ref_all(:, wl), 'omitnan');
end

%  ==============================
%%      All TD-NIRS files
%  ==============================

files = dir('DPF_MmAB*.DAT');
nFiles = numel(files);

if nFiles == 0
    error('Nessun file DPF_MmAB*.DAT trovato nella cartella corrente.');
end

% Results
DPF_all_reps  = zeros(nFiles, 5, 4);   % file x rip x lambda
DPF_mean_wl   = zeros(nFiles, 4);      % file x lambda
DPF_std_wl    = zeros(nFiles, 4);      % file x lambda

for f = 1:nFiles
    Settings.inputFile = files(f).name;

    % Reading the curves for the current file
    curves = TDNIRS_dataLoad(Settings);   % [4096 x 20]
    curves = double(curves);

    % Reshape in [4096 x 5 x 4] - (5 rip x 4 lambda, checked using TDNIRS_dataLoad)
    curves3D = reshape(curves, Settings.nChannels, 5, 4);

    for wl = 1:4
        for rep = 1:5
            y = curves3D(:, rep, wl);

            % Baseline correction
            y = y - min(y);
            y(y < 0) = 0;

            % Avoid dividing by 0
            if sum(y) == 0
                DPF_all_reps(f, rep, wl) = NaN;
                continue;
            end

            % Mean time of flight sample
            mean_t = sum(t .* y) / sum(y);

            % DPF with IRF reference
            DPF_all_reps(f, rep, wl) = c0 * (mean_t - mean_t_ref_wl(wl)) / (n * rho);
        end

        % Mean and std over the 5 repetitions
        DPF_mean_wl(f, wl) = mean(DPF_all_reps(f, :, wl), 'omitnan');
        DPF_std_wl(f, wl)  = std(DPF_all_reps(f, :, wl),  'omitnan');
    end

    fprintf('File %s elaborato.\n', files(f).name);
end

% Mean and std across all placements for each wavelength
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

%   ==================================
%%      Saving and Visualizing Results
%   ==================================

T = table( ...
    string({files.name})', ...
    DPF_mean_wl(:,1), DPF_std_wl(:,1), ...
    DPF_mean_wl(:,2), DPF_std_wl(:,2), ...
    DPF_mean_wl(:,3), DPF_std_wl(:,3), ...
    DPF_mean_wl(:,4), DPF_std_wl(:,4), ...
    'VariableNames', { ...
    'File', ...
    'DPF_670_mean','DPF_670_std', ...
    'DPF_730_mean','DPF_730_std', ...
    'DPF_780_mean','DPF_780_std', ...
    'DPF_830_mean','DPF_830_std'});

disp(T);
disp(T_allplacements);

fprintf('\nEstimated DPF at 760 nm (local linear interpolation 730-780): %.6f\n', DPF_760);
fprintf('Estimated DPF at 850 nm (linear extrapolation 780-830): %.6f\n', DPF_850);

% % Save CSV
% writetable(T, 'DPF_results_allfiles.csv');
% disp('Results saved in DPF_results_allfiles.csv');

%% DPF heatmap per file and wavelength
figure;
h = heatmap(string(lambda) + " nm", string({files.name})', DPF_mean_wl, ...
    'CellLabelFormat', '%.3f');
gcf.Color = 'w';
h.FontColor = 'k';
h.Title = 'Mean DPF heatmap per file and wavelength';
h.XLabel = 'Wavelength';
h.YLabel = 'File';
h.Colormap = jet(256);
h.ColorLimits = [min(DPF_mean_wl, [], 'all', 'omitnan'), max(DPF_mean_wl, [], 'all', 'omitnan')];
