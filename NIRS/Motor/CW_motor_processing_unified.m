%% HOMER3 CW-fNIRS pipeline using native Homer3 Block Average
% This script:
% 1) Loads a SNIRF file using Homer3
% 2) Converts intensity to OD
% 3) Optionally performs motion artifact detection
% 4) Optionally performs wavelet correction
% 5) Optionally performs spline correction
% 6) Optionally performs combined wavelet + spline correction
% 7) Optionally performs bandpass filtering
% 7) Converts OD to concentration (HbO / HbR / HbT)
% 8) Computes native Homer3 Block Average using hmrR_BlockAvg
% 9) Plots single-channel and all-channel Block Average
% 10) Builds HRF heatmap and HRF 3D surface
% 11) Plots probe geometry and an HRF metric on the probe

close all; clc

%% ================= PATHS =================
HOMER3_DIR = '/home/giulia/MATLAB/Toolboxes/Homer3';
addpath(genpath(HOMER3_DIR));

%% ================= INPUT =================
snirfFile = '/home/giulia/Scrivania/NIRS multimodal/Motor/Motor/CW - Motor/2025-04-01_001.snirf';

keepStims = ["1","2"];
% Only these stimulus labels will be kept if they exist in the file.

%% ================= USER OPTIONS =================
exampleCh   = 1;   % Example S-D pair index used for single-channel plots
hbTypeIdx   = 1;   % 1 = HbO, 2 = HbR, 3 = HbT
condIdxPlot = 1;   % Condition used for all-channel BlockAvg, heatmap, 3D plots

plotAllBlockAvgChannels = true;
channelsPerFigure = 9;

showProbe3D            = true;
showProbe3D_HRF        = true;
probeMetricHbTypeIdx   = 1;      % 1 = HbO, 2 = HbR, 3 = HbT
probeMetricTimeWindow  = [4 8];  % Time window (s) for probe HRF metric
probeUseWavelengthIdx  = 1;      % Wavelength index used to map channels on probe
showChannelLinks       = true;

%% ================= DENOISING OPTIONS =================
useMotionArtifactDetect    = true;
useMotionArtifactByChannel = false;

if ~exist('motionCorrectionMethod', 'var') || isempty(motionCorrectionMethod)
    motionCorrectionMethod = "wavelet";
end
motionCorrectionMethod = string(lower(motionCorrectionMethod));

% Options:
% "none"
% "wavelet"
% "spline"
% "combined"

validMotionCorrectionMethods = ["none", "wavelet", "spline", "combined"];
assert(any(motionCorrectionMethod == validMotionCorrectionMethods), ...
    'Unknown motionCorrectionMethod: %s', motionCorrectionMethod);

if any(motionCorrectionMethod == ["spline", "combined"]) && ~useMotionArtifactByChannel
    fprintf(['Enabling channel-wise motion artifact detection because ' ...
        'method "%s" requires tIncCh.\n'], motionCorrectionMethod);
    useMotionArtifactByChannel = true;
end

useBandpassFilter = true;

% ---- Motion artifact detection parameters ----
tMotion     = 0.5;
tMask       = 1.0;
STDEVthresh = 8;
AMPthresh   = 1;

% ---- Wavelet parameters ----
iqr_wavelet    = 1.0;  %"Wavelet-based motion artifact removal for functional near-infrared spectroscopy"
turnon_wavelet = 1;

% ---- Spline parameters ----
pSpline = 0.01; %"How to detect and reduce movement artifacts in near-infrared imaging using moving standard deviation and spline interpolation"

% ---- Bandpass parameters ----
hpf = 0.01;
lpf = 0.2;

% ---- MBLL parameters ----
dpf = [5.801187 5.119572]; %from interpolation of TD DPF values (used for reference)

%% ================= BLOCK AVERAGE PARAMETERS =================
trangeBA = [-2 20];

%% ================= SAVE =================
outMat = sprintf('CW_results_%s.mat', motionCorrectionMethod);
%% ================= LOAD SNIRF =================
assert(exist(snirfFile,'file')==2, 'SNIRF not found: %s', snirfFile);
% Check that the file exists before trying to load it.
% Input:
%   snirfFile -> full path to the .snirf file
% Output:
%   none
% If the file is missing, execution stops here.

snirf   = SnirfClass(snirfFile);
% Homer3 class constructor.
% Input:
%   snirfFile -> path to a SNIRF file
% Output:
%   snirf -> a SNIRF container object with fields such as:
%            - data
%            - probe
%            - stim
%            - aux (if present)
%
% This is the main object that represents the recording.

dataObj = snirf.data;
% Extract the SNIRF.data container from the file.
% Input:
%   snirf -> SNIRF container
% Output:
%   dataObj -> Homer3 data container
%
% This contains the actual recorded time series.

probe   = snirf.probe;
% Extract the probe definition.
% Input:
%   snirf
% Output:
%   probe -> Homer3 ProbeClass object
%
% The probe stores:
%   - source positions
%   - detector positions
%   - wavelengths
%   - landmarks
%   - coordinate system


% Extract the first data object robustly
% Why needed:
% Homer3 object storage can vary depending on file/version.
d0 = get_first_data_element(dataObj);

% Extract time and numeric data
t  = get_time_from_data(d0);

% Custom helper.
% Input:
%   d0 -> one Homer3 data object
% Output:
%   Y0 -> numeric matrix [time x measurements]
%
% The helper tries:
%   - GetDataTimeSeries()
%   - .dataTimeSeries
%   - GetDataTimeSeries('reshape')
Y0 = get_data_matrix(d0);

nMeas = size(Y0, 2);

% Extract measurement list
% Typical columns in raw OD/intensity measurement list:
%   col 1 -> source index
%   col 2 -> detector index
%   col 3 -> data type / placeholder
%   col 4 -> wavelength index
ml0 = get_meas_list_safe(d0);



%% ========= RAW SIGNAL QUALITY CHECK - question 1 ==========
meanRaw = mean(Y0,1,'omitnan');
stdRaw  = std(Y0,0,1,'omitnan');
cvRaw   = stdRaw ./ meanRaw;

badIntensity = meanRaw < prctile(meanRaw,5);
badCV = cvRaw > 0.05;   
badChannels  = badIntensity | badCV;

fprintf('Potential bad measurements: %d / %d\n', sum(badChannels), numel(badChannels));
disp(find(badChannels))

figure('Name','Raw intensity QC','Color','w');
subplot(2,1,1)
bar(meanRaw)
xlabel('Measurement'); ylabel('Mean intensity')
title('Mean raw intensity per measurement')
grid on

subplot(2,1,2)
bar(cvRaw)
xlabel('Measurement'); ylabel('Coefficient of variation')
title('Relative variability per measurement')
grid on

figure('Name','All channels with BAD highlighted','Color','w'); 
hold on

goodIdx = find(~badChannels);
badIdx  = find(badChannels);

% ---- GOOD channels (grigio chiaro) ----
for i = goodIdx
    plot(t, Y0(:,i), 'Color', [0.7 0.7 0.7])
end

% ---- BAD channels (rosso, più spessi) ----
for i = badIdx
    plot(t, Y0(:,i), 'r', 'LineWidth', 1.5)
end

xlabel('Time (s)')
ylabel('Intensity')
title('All channels (gray) vs BAD channels (red)')
grid on

legend({'Good','Bad'})


%%
% Build active measurement lists in Homer3 format:
% [source detector active wavelength_or_type]
if size(ml0,2) >= 4
    mlActMan  = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,4)] };
    mlActAuto = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,4)] };
elseif size(ml0,2) >= 3
    mlActMan  = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,3)] };
    mlActAuto = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,3)] };
else
    error('Unexpected measurement list format in raw data.');
end

fprintf('Loaded: %s\n', snirfFile);
fprintf('Time points: %d | nMeas: %d\n', numel(t), nMeas);
fprintf('Wavelengths: %s\n', mat2str(get_probe_wavelengths(probe)));
fprintf('MeasList size: %d x %d\n', size(ml0,1), size(ml0,2));

%% ================= KEEP STIMS =================
hasStim = false;
try
    hasStim = ~isempty(snirf.stim);
catch
end

if hasStim
    stimNamesAll = get_stim_names(snirf.stim);
    % Custom helper.
    % Input:
    %   snirf.stim -> SNIRF stimulus container
    % Output:
    %   stimNamesAll -> string array with all condition labels
    fprintf('Stim names found in file:\n');
    disp(stimNamesAll)

    % Keep requested stimuli AND force the order given in keepStims
    stimNamesClean = lower(strtrim(stimNamesAll));
    keepStimsClean = lower(strtrim(keepStims));

    [tf, ord] = ismember(keepStimsClean, stimNamesClean);

    if any(tf)
        ord = ord(tf);
        snirf.stim = snirf.stim(ord);
        fprintf('Kept stims in requested order: %s\n', strjoin(get_stim_names(snirf.stim), ', '));
    else
        warning('Requested stimulus labels were not found. Using all available stimuli.');
        fprintf('Available stimuli: %s\n', strjoin(stimNamesAll, ', '));
    end
else
    warning('No stimuli found in SNIRF.');
end

condNames = get_stim_names(snirf.stim);


fprintf('Final number of conditions: %d\n', numel(condNames));
disp(condNames)

if isempty(condNames)
    error('No valid conditions available after stimulus selection.');
end

if condIdxPlot > numel(condNames)
    warning('condIdxPlot=%d is out of range. Resetting to 1.', condIdxPlot);
    condIdxPlot = 1;
end

%% ================= PIPELINE =================
% -------- 1) Intensity -> Optical Density --------
dodRawObj = hmrR_Intensity2OD(dataObj);

% -------- 2) Motion artifact detection (optional) --------
tInc   = [];
tIncCh = [];
tIncGlobalFromChannel = [];
tIncMan = {ones(length(t),1)};

if useMotionArtifactDetect
    try
        tInc = hmrR_MotionArtifact(dodRawObj, probe, mlActMan, mlActAuto, tIncMan, ...
                                   tMotion, tMask, STDEVthresh, AMPthresh);
    catch ME
        warning('Global motion artifact detection failed: %s', ME.message);
        tInc = [];
    end
end


if useMotionArtifactByChannel
    try
        [tIncGlobalFromChannel, tIncCh] = hmrR_MotionArtifactByChannel( ...
            dodRawObj, probe, mlActMan, mlActAuto, tIncMan, ...
            tMotion, tMask, STDEVthresh, AMPthresh);

        if isempty(tInc) && ~isempty(tIncGlobalFromChannel)
            tInc = tIncGlobalFromChannel;
        end
    catch ME
        warning('Channel-wise motion artifact detection failed: %s', ME.message);
        tIncGlobalFromChannel = [];
        tIncCh = [];
    end
end

%% ================= MOTION ARTIFACT SUMMARY - question 2
if ~isempty(tInc)
    if iscell(tInc)
        tIncNum = tInc{1};
    else
        tIncNum = tInc;
    end

    rejectedPercent = 100 * mean(tIncNum == 0);
    fprintf('Motion artifact rejected samples: %.2f %%\n', rejectedPercent);

    figure('Name','Motion artifact mask','Color','w');
    plot(t, tIncNum, 'k', 'LineWidth', 1.2)
    ylim([-0.1 1.1])
    xlabel('Time (s)')
    ylabel('tInc')
    title(sprintf('Motion artifact inclusion mask | rejected %.2f %%', rejectedPercent))
    grid on
end
%%

% -------- 3) Motion correction: choose one method --------

dodCorrObj = dodRawObj;
dodWaveObj = dodRawObj;
dodSplineObj = dodRawObj;

switch lower(motionCorrectionMethod)

    case "none"
        fprintf('No motion correction applied.\n');
        dodCorrObj = dodRawObj;

    case "wavelet"
        fprintf('Applying wavelet motion correction only.\n');

        dodWaveObj = hmrR_MotionCorrectWavelet(dodRawObj, mlActMan, mlActAuto, ...
                                               iqr_wavelet, turnon_wavelet);

        dodCorrObj = dodWaveObj;

    case "spline"
        fprintf('Applying spline motion correction only.\n');

        try
            tIncCh_use = build_spline_tIncCh(dodRawObj, tInc, tIncCh);
            dodSplineObj = hmrR_MotionCorrectSpline( ...
                dodRawObj, mlActAuto, tIncCh_use, pSpline, 1);

            dodCorrObj = dodSplineObj;

        catch ME
            warning('Spline correction failed: %s. Using raw OD instead.', ME.message);
            dodCorrObj = dodRawObj;
        end

    case "combined"
        fprintf('Applying combined wavelet + spline motion correction.\n');

        try
            dodWaveObj = hmrR_MotionCorrectWavelet(dodRawObj, mlActMan, mlActAuto, ...
                                                   iqr_wavelet, turnon_wavelet);
            tIncCh_use = build_spline_tIncCh(dodRawObj, tInc, tIncCh);
            dodSplineObj = hmrR_MotionCorrectSpline( ...
                dodWaveObj, mlActAuto, tIncCh_use, pSpline, 1);

            dodCorrObj = dodSplineObj;

        catch ME
            warning('Combined correction failed: %s. Falling back to wavelet only.', ME.message);
            dodCorrObj = dodWaveObj;
        end

    otherwise
        error('Unknown motionCorrectionMethod: %s', motionCorrectionMethod);
end
%

% -------- 4) Bandpass filtering (optional) --------
if useBandpassFilter
    dodFiltObj = hmrR_BandpassFilt(dodCorrObj, hpf, lpf);
else
    dodFiltObj = dodCorrObj;
end

% -------- 5) OD -> concentration --------
wl = get_probe_wavelengths(probe);
assert(numel(dpf)==numel(wl), ...
    'dpf length (%d) must match number of wavelengths (%d).', numel(dpf), numel(wl));

dcObj = hmrR_OD2Conc(dodFiltObj, probe, dpf); %dpf same for all channels and wavelengths

%% ================= EXTRACT OD DATA =================
dodRaw    = get_data_matrix(get_first_data_element(dodRawObj));
dodWave   = get_data_matrix(get_first_data_element(dodWaveObj));
dodSpline = get_data_matrix(get_first_data_element(dodSplineObj));
dodFilt   = get_data_matrix(get_first_data_element(dodFiltObj));
dodCorr = get_data_matrix(get_first_data_element(dodCorrObj));
ml        = get_meas_list_safe(get_first_data_element(dodRawObj));

nCh = size(dodRaw,2);
assert(exampleCh >= 1 && exampleCh <= max(1, estimate_num_sd_pairs(ml)), ...
    'exampleCh out of range.');

fprintf('OD raw size      : %d x %d\n', size(dodRaw,1), size(dodRaw,2));
fprintf('OD wavelet size  : %d x %d\n', size(dodWave,1), size(dodWave,2));
fprintf('OD spline size   : %d x %d\n', size(dodSpline,1), size(dodSpline,2));
fprintf('OD filtered size : %d x %d\n', size(dodFilt,1), size(dodFilt,2));

%% ================= OD PLOT: CHOSEN CHANNEL =================
figure('Name', sprintf('OD example measurement %d', exampleCh), 'Color', 'w');

plot(t, dodRaw(:,exampleCh),  'k-', 'LineWidth', 1); hold on;
plot(t, dodCorr(:,exampleCh), 'b-', 'LineWidth', 1.1);
plot(t, dodFilt(:,exampleCh), 'r-', 'LineWidth', 1.3);

grid on
xlabel('Time (s)');
ylabel('OD');
title(sprintf('OD example measurement | %s | method: %s', ...
    make_channel_title(ml, exampleCh), motionCorrectionMethod), ...
    'Interpreter', 'none');

legend({'Raw OD','After motion correction','After bandpass'}, 'Location', 'best');
%% ================= OD PLOT: RAW VS FINAL =================
figure('Name', sprintf('OD raw vs filtered - measurement %d', exampleCh), 'Color', 'w');
plot(t, dodRaw(:,exampleCh),  'k-', 'LineWidth', 1); hold on;
plot(t, dodFilt(:,exampleCh), 'r-', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('OD');
title(sprintf('OD raw vs final filtered | %s', make_channel_title(ml, exampleCh)), 'Interpreter', 'none');
legend({'Raw OD','Final filtered OD'}, 'Location', 'best');

%% ================= OD HEATMAPS =================
figure('Name','Raw OD heatmap','Color','w');
imagesc(t, 1:nCh, dodRaw');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Raw OD - all measurements');
colorbar

figure('Name','Wavelet-stage OD heatmap','Color','w');
imagesc(t, 1:nCh, dodWave');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Wavelet-stage OD - all measurements');
colorbar

figure('Name','Spline-stage OD heatmap','Color','w');
imagesc(t, 1:nCh, dodSpline');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Spline-stage OD - all measurements');
colorbar

figure('Name','Final filtered OD heatmap','Color','w');
imagesc(t, 1:nCh, dodFilt');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Final filtered OD - all measurements');
colorbar

%% ================= NATIVE HOMER3 BLOCK AVERAGE =================
% Each condition is processed separately using the native Homer3 stim container.
% The output is a Homer3 data object. We then extract:
% - numeric matrix
% - time vector
% - measurement list
% and reconstruct HbO/HbR/HbT channel-by-channel.

HRF_by_condition  = cell(numel(condNames),1);
tHRF_by_condition = cell(numel(condNames),1);
mlBA_by_condition = cell(numel(condNames),1);

figure('Name','BlockAvg single channel HbO/HbR','Color','w');
nCond = numel(condNames);

for c = 1:numel(condNames)

    stim_single = snirf.stim(c);

    [dataAvgBA, dataStdBA, nTrialsBA, dataSum2BA, yTrialsBA] = ...
        hmrR_BlockAvg(dcObj, stim_single, trangeBA);

    yBA = get_data_matrix(get_first_data_element(dataAvgBA));
    tHRF_BA = get_time_from_data(get_first_data_element(dataAvgBA));
    mlBA = get_meas_list_safe(get_first_data_element(dataAvgBA));

    HRF_by_condition{c}  = yBA;
    tHRF_by_condition{c} = tHRF_BA(:);
    mlBA_by_condition{c} = mlBA;

    subplot(nCond,1,c); hold on; grid on

    % Reconstruct HbO/HbR/HbT for the chosen S-D pair index
    idxHbO = find_hb_channels(mlBA, exampleCh, 'HbO');
    idxHbR = find_hb_channels(mlBA, exampleCh, 'HbR');
    idxHbT = find_hb_channels(mlBA, exampleCh, 'HbT');

    if ~isempty(idxHbO)
        plot(tHRF_BA, yBA(:,idxHbO(1)), 'r-', 'LineWidth', 2);
    end
    if ~isempty(idxHbR)
        plot(tHRF_BA, yBA(:,idxHbR(1)), 'b-', 'LineWidth', 2);
    end
    if ~isempty(idxHbT)
        plot(tHRF_BA, yBA(:,idxHbT(1)), 'k--', 'LineWidth', 1.5);
        legend({'HbO','HbR','HbT'}, 'Location', 'best');
    else
        legend({'HbO','HbR'}, 'Location', 'best');
    end

    xlabel('Time (s)');
    ylabel('\DeltaHb (\muM)');
    title(sprintf('BlockAvg | %s | channel %d', condNames(c), exampleCh), ...
        'Interpreter','none');
end


%% ================= HRF METRICS PER TASK - question 3 ========
metricWindow = [4 8];

Condition = strings(0,1);
Channel   = [];
HbO_mean  = [];
HbO_peak  = [];
HbR_mean  = [];
HbR_min   = [];

for c = 1:numel(condNames)

    yBA = HRF_by_condition{c};
    tH  = tHRF_by_condition{c};
    mlB = mlBA_by_condition{c};

    if isempty(yBA)
        continue
    end

    nPairs = estimate_num_sd_pairs(mlB);
    mask = tH >= metricWindow(1) & tH <= metricWindow(2);

    for ch = 1:nPairs

        idxHbO = find_hb_channels(mlB, ch, 'HbO');
        idxHbR = find_hb_channels(mlB, ch, 'HbR');

        Condition(end+1,1) = condNames(c);
        Channel(end+1,1)   = ch;

        if ~isempty(idxHbO)
            hbo = yBA(:,idxHbO(1));
            HbO_mean(end+1,1) = mean(hbo(mask), 'omitnan');
            HbO_peak(end+1,1) = max(hbo(mask), [], 'omitnan');
        else
            HbO_mean(end+1,1) = NaN;
            HbO_peak(end+1,1) = NaN;
        end

        if ~isempty(idxHbR)
            hbr = yBA(:,idxHbR(1));
            HbR_mean(end+1,1) = mean(hbr(mask), 'omitnan');
            HbR_min(end+1,1)  = min(hbr(mask), [], 'omitnan');
        else
            HbR_mean(end+1,1) = NaN;
            HbR_min(end+1,1)  = NaN;
        end
    end
end

Results = table(Condition, Channel, HbO_mean, HbO_peak, HbR_mean, HbR_min);

Results.Properties.VariableNames = { ...
    'Condition', 'Channel', ...
    'HbO_mean_4_8', 'HbO_peak_4_8', ...
    'HbR_mean_4_8', 'HbR_min_4_8'};

% REMOVE BAD CHANNELS 
badMeasIdx = find(badChannels);
badPairs = unique(ml0(badMeasIdx,1:2), 'rows');

allPairs = unique(ml0(:,1:2), 'rows', 'stable');
badPairIdx = find(ismember(allPairs, badPairs, 'rows'));

Results = Results(~ismember(Results.Channel, badPairIdx), :);

disp(Results)
writetable(Results, 'CW_HRF_metrics_by_task.csv');


%% ================= HEMISPHERIC ANALYSIS: LEFT vs RIGHT =================
% Task 1 = left hand grip  -> expected right hemisphere activation
% Task 2 = right hand grip -> expected left hemisphere activation

srcPos = get_probe_srcPos(probe);
detPos = get_probe_detPos(probe);

pairMap = unique(ml0(:,1:2), 'rows', 'stable');
idxPairs = 1:size(pairMap,1);

[midPos, srcDetPairs] = compute_pair_midpoints(pairMap, srcPos, detPos, idxPairs);

% IMPORTANT:
% Here we assume X < 0 = left hemisphere, X > 0 = right hemisphere.
% If your probe plot shows the opposite, switch the signs below.
hemiLabel = strings(size(midPos,1),1);
hemiLabel(midPos(:,1) < 0) = "Left";
hemiLabel(midPos(:,1) > 0) = "Right";
hemiLabel(midPos(:,1) == 0) = "Midline";

% Add hemisphere info to Results
Results.Hemisphere = strings(height(Results),1);

for r = 1:height(Results)
    ch = Results.Channel(r);
    if ch <= numel(hemiLabel)
        Results.Hemisphere(r) = hemiLabel(ch);
    else
        Results.Hemisphere(r) = "Unknown";
    end
end

% Remove midline/unknown if present
ResultsHemi = Results(Results.Hemisphere == "Left" | Results.Hemisphere == "Right", :);

% Separate task 1 and task 2
R_task1 = ResultsHemi(ResultsHemi.Condition == "1", :);  % left hand
R_task2 = ResultsHemi(ResultsHemi.Condition == "2", :);  % right hand

% HbO mean per hemisphere
fprintf('\n===== HEMISPHERIC ACTIVATION ANALYSIS =====\n');

fprintf('\nTask 1 = LEFT hand grip\n');
mean_L_task1 = mean(R_task1.HbO_mean_4_8(R_task1.Hemisphere=="Left"), 'omitnan');
mean_R_task1 = mean(R_task1.HbO_mean_4_8(R_task1.Hemisphere=="Right"), 'omitnan');

fprintf('Left hemisphere HbO : %.4e\n', mean_L_task1);
fprintf('Right hemisphere HbO: %.4e\n', mean_R_task1);

fprintf('\nTask 2 = RIGHT hand grip\n');
mean_L_task2 = mean(R_task2.HbO_mean_4_8(R_task2.Hemisphere=="Left"), 'omitnan');
mean_R_task2 = mean(R_task2.HbO_mean_4_8(R_task2.Hemisphere=="Right"), 'omitnan');

fprintf('Left hemisphere HbO : %.4e\n', mean_L_task2);
fprintf('Right hemisphere HbO: %.4e\n', mean_R_task2);

% Statistical comparison: Left vs Right within each task
HbO_L_task1 = R_task1.HbO_mean_4_8(R_task1.Hemisphere=="Left");
HbO_R_task1 = R_task1.HbO_mean_4_8(R_task1.Hemisphere=="Right");

HbO_L_task2 = R_task2.HbO_mean_4_8(R_task2.Hemisphere=="Left");
HbO_R_task2 = R_task2.HbO_mean_4_8(R_task2.Hemisphere=="Right");

[~, p_task1_hemi] = ttest2(HbO_L_task1, HbO_R_task1);
[~, p_task2_hemi] = ttest2(HbO_L_task2, HbO_R_task2);

fprintf('\nTask 1 Left vs Right hemisphere HbO: p = %.4f\n', p_task1_hemi);
fprintf('Task 2 Left vs Right hemisphere HbO: p = %.4f\n', p_task2_hemi);

% Plot hemispheric comparison
figure('Name','Hemispheric comparison HbO','Color','w');

subplot(1,2,1)
boxplot([HbO_L_task1; HbO_R_task1], ...
        [repmat("Left",numel(HbO_L_task1),1); repmat("Right",numel(HbO_R_task1),1)])
ylabel('Mean HbO 4-8 s')
title(sprintf('Task 1 left hand | p = %.4f', p_task1_hemi))
grid on

subplot(1,2,2)
boxplot([HbO_L_task2; HbO_R_task2], ...
        [repmat("Left",numel(HbO_L_task2),1); repmat("Right",numel(HbO_R_task2),1)])
ylabel('Mean HbO 4-8 s')
title(sprintf('Task 2 right hand | p = %.4f', p_task2_hemi))
grid on


%% ===== Hemispheric HbO comparison - question=====

vals = {HbO_L_task1, HbO_R_task1, HbO_L_task2, HbO_R_task2};

means = [ ...
    mean(HbO_L_task1,'omitnan'), mean(HbO_R_task1,'omitnan'); ...
    mean(HbO_L_task2,'omitnan'), mean(HbO_R_task2,'omitnan') ...
];

sems = [ ...
    std(HbO_L_task1,'omitnan')/sqrt(numel(HbO_L_task1)), std(HbO_R_task1,'omitnan')/sqrt(numel(HbO_R_task1)); ...
    std(HbO_L_task2,'omitnan')/sqrt(numel(HbO_L_task2)), std(HbO_R_task2,'omitnan')/sqrt(numel(HbO_R_task2)) ...
];

figure('Name','Hemispheric activation HbO - nice plot','Color','w');
hold on

b = bar(means);
b(1).FaceColor = [0.35 0.55 0.90];   % Left hemisphere
b(2).FaceColor = [0.90 0.45 0.35];   % Right hemisphere

% Error bars
ngroups = size(means,1);
nbars = size(means,2);
x = nan(nbars, ngroups);

for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

errorbar(x', means, sems, 'k', ...
    'linestyle','none', 'LineWidth',1.2, 'CapSize',10);


xticks([1 2])
xticklabels({'Task 1: left hand','Task 2: right hand'})
ylabel('Mean HbO 4-8 s')
title('Hemispheric activation during left and right hand grip')
legend({'Left hemisphere','Right hemisphere'}, 'Location','best')
grid on
box off

% Add p-values above bars
yl = ylim;
yText1 = max(means(1,:) + sems(1,:)) + 0.08*range(yl);
yText2 = max(means(2,:) + sems(2,:)) + 0.08*range(yl);

text(1, yText1, sprintf('p = %.4f', p_task1_hemi), ...
    'HorizontalAlignment','center', 'FontWeight','bold')

text(2, yText2, sprintf('p = %.4f', p_task2_hemi), ...
    'HorizontalAlignment','center', 'FontWeight','bold')
%% ================= SNR CALCULATION =================
baselineWindow = [-2 0];   % prima dello stimolo
signalWindow   = [4 8];    % risposta HRF

Condition = strings(0,1);
Channel   = [];
SNR_HbO   = [];

for c = 1:numel(condNames)

    yBA = HRF_by_condition{c};
    tH  = tHRF_by_condition{c};
    mlB = mlBA_by_condition{c};

    if isempty(yBA)
        continue
    end

    nPairs = estimate_num_sd_pairs(mlB);

    baseMask = tH >= baselineWindow(1) & tH <= baselineWindow(2);
    sigMask  = tH >= signalWindow(1)   & tH <= signalWindow(2);

    for ch = 1:nPairs

        idxHbO = find_hb_channels(mlB, ch, 'HbO');

        if isempty(idxHbO)
            continue
        end

        hbo = yBA(:, idxHbO(1));

        baselineStd = std(hbo(baseMask), 'omitnan');
        signalAmp   = mean(hbo(sigMask), 'omitnan');

        Condition(end+1,1) = condNames(c);
        Channel(end+1,1)   = ch;
        SNR_HbO(end+1,1)   = signalAmp / baselineStd;
    end
end

SNR_table = table(Condition, Channel, SNR_HbO);

SNR_table = SNR_table(~ismember(SNR_table.Channel, badPairIdx), :);
disp(SNR_table)
%% ================= TASK COMPARISON: LEFT VS RIGHT GRIP - question 4 =================
% NOTE:
% Single-subject analysis → no inferential statistics (no t-test).
% Channels are not independent observations.

if numel(condNames) >= 2

    cond1 = condNames(1);
    cond2 = condNames(2);

    R1 = Results(Results.Condition == cond1, :);
    R2 = Results(Results.Condition == cond2, :);

    % Align channels
    [commonCh, idx1, idx2] = intersect(R1.Channel, R2.Channel);

    HbO_1 = R1.HbO_mean_4_8(idx1);
    HbO_2 = R2.HbO_mean_4_8(idx2);

    HbR_1 = R1.HbR_mean_4_8(idx1);
    HbR_2 = R2.HbR_mean_4_8(idx2);

    % ---- DESCRIPTIVE STATS ----
    mean_HbO_1 = mean(HbO_1, 'omitnan');
    mean_HbO_2 = mean(HbO_2, 'omitnan');

    mean_HbR_1 = mean(HbR_1, 'omitnan');
    mean_HbR_2 = mean(HbR_2, 'omitnan');

    fprintf('\n===== TASK COMPARISON (DESCRIPTIVE) =====\n');
    fprintf('%s HbO mean: %.4e\n', cond1, mean_HbO_1);
    fprintf('%s HbO mean: %.4e\n', cond2, mean_HbO_2);

    fprintf('%s HbR mean: %.4e\n', cond1, mean_HbR_1);
    fprintf('%s HbR mean: %.4e\n', cond2, mean_HbR_2);

    % ---- PLOT ----
    figure('Name','Grip task comparison HbO (descriptive)','Color','w');
    boxplot([HbO_1 HbO_2], 'Labels', cellstr([cond1 cond2]))
    ylabel('Mean HbO 4-8 s')
    title(sprintf('HbO comparison: %s vs %s (descriptive)', cond1, cond2), ...
        'Interpreter','none')
    grid on

    figure('Name','Grip task comparison HbR (descriptive)','Color','w');
    boxplot([HbR_1 HbR_2], 'Labels', cellstr([cond1 cond2]))
    ylabel('Mean HbR 4-8 s')
    title(sprintf('HbR comparison: %s vs %s (descriptive)', cond1, cond2), ...
        'Interpreter','none')
    grid on
end

%% ================= BLOCKAVG PLOT FOR ALL CHANNELS =================
if plotAllBlockAvgChannels
    yBA_plot  = HRF_by_condition{condIdxPlot};
    tHRF_plot = tHRF_by_condition{condIdxPlot};
    mlBA_plot = mlBA_by_condition{condIdxPlot};

    if isempty(yBA_plot)
        warning('Selected condition has no BlockAvg output.');
    else
        nPairs = estimate_num_sd_pairs(mlBA_plot);

        nFig = ceil(nPairs / channelsPerFigure);
        nRow = ceil(sqrt(channelsPerFigure));
        nCol = ceil(channelsPerFigure / nRow);

        for f = 1:nFig
            figure('Name', sprintf('BlockAvg channels %d', f), 'Color', 'w');

            chStart = (f-1)*channelsPerFigure + 1;
            chEnd   = min(f*channelsPerFigure, nPairs);

            for ch = chStart:chEnd
                subplot(nRow, nCol, ch-chStart+1)
                hold on; grid on

                idxHbO = find_hb_channels(mlBA_plot, ch, 'HbO');
                idxHbR = find_hb_channels(mlBA_plot, ch, 'HbR');
                idxHbT = find_hb_channels(mlBA_plot, ch, 'HbT');

                if ~isempty(idxHbO)
                    plot(tHRF_plot, yBA_plot(:,idxHbO(1)), 'r-', 'LineWidth', 1.5);
                end
                if ~isempty(idxHbR)
                    plot(tHRF_plot, yBA_plot(:,idxHbR(1)), 'b-', 'LineWidth', 1.5);
                end
                if ~isempty(idxHbT)
                    plot(tHRF_plot, yBA_plot(:,idxHbT(1)), 'k--', 'LineWidth', 1.0);
                end

                xlabel('Time (s)');
                ylabel('\DeltaHb');
                title(sprintf('Ch %d', ch));
            end

            sgtitle(sprintf('BlockAvg per channel | %s', condNames(condIdxPlot)), ...
                'Interpreter', 'none');
        end
    end
end

%% ================= HRF HEATMAP + HRF 3D =================
yBA_plot  = HRF_by_condition{condIdxPlot};
tHRF_plot = tHRF_by_condition{condIdxPlot};
mlBA_plot = mlBA_by_condition{condIdxPlot};

if isempty(yBA_plot)
    warning('HRF for selected condition is empty. Skipping HRF heatmap and 3D plot.');
else
    hbName = hb_idx_to_name(hbTypeIdx);

    nPairs = estimate_num_sd_pairs(mlBA_plot);
    HRF_all = nan(length(tHRF_plot), nPairs);

    for ch = 1:nPairs
        idxHb = find_hb_channels(mlBA_plot, ch, hbName);
        if ~isempty(idxHb)
            HRF_all(:,ch) = yBA_plot(:,idxHb(1));
        end
    end

    figure('Name','HRF heatmap','Color','w');
    imagesc(tHRF_plot, 1:nPairs, HRF_all');
    axis xy
    xlabel('Time (s)');
    ylabel('Channel');
    title(sprintf('HRF heatmap | %s | %s', condNames(condIdxPlot), hbName), ...
        'Interpreter','none');
    colorbar

    [T_hrf, Ch_hrf] = meshgrid(tHRF_plot, 1:nPairs);

    figure('Name','HRF 3D surface','Color','w');
    surf(T_hrf, Ch_hrf, HRF_all', 'EdgeColor', 'none');
    xlabel('Time (s)');
    ylabel('Channel');
    zlabel('\DeltaHb (\muM)');
    title(sprintf('HRF 3D | %s | %s', condNames(condIdxPlot), hbName), ...
        'Interpreter','none');
    colorbar
    view(45,30)
    grid on
end

%% ================= PROBE GEOMETRY 2D/3D =================
if showProbe3D
    srcPos  = get_probe_srcPos(probe);
    detPos  = get_probe_detPos(probe);
    landPos = get_probe_landmarkPos(probe);

    if isempty(srcPos) || isempty(detPos)
        warning('Probe plot skipped: source/detector coordinates not available.');
    else
        is3D = (size(srcPos,2) >= 3) && (size(detPos,2) >= 3);

        if is3D
            figure('Name','Probe geometry 3D','Color','w'); hold on

            if ~isempty(landPos) && size(landPos,2) >= 3
                scatter3(landPos(:,1), landPos(:,2), landPos(:,3), 30, [0.7 0.7 0.7], 'filled');
            end

            scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 90, 'r', 'filled');
            scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 90, 'b', 'filled');

            for i = 1:size(srcPos,1)
                text(srcPos(i,1), srcPos(i,2), srcPos(i,3), sprintf(' S%d', i), ...
                    'Color', 'r', 'FontSize', 9);
            end

            for i = 1:size(detPos,1)
                text(detPos(i,1), detPos(i,2), detPos(i,3), sprintf(' D%d', i), ...
                    'Color', 'b', 'FontSize', 9);
            end

            if showChannelLinks && ~isempty(ml)
                for ch = 1:size(ml,1)
                    sIdx = ml(ch,1);
                    dIdx = ml(ch,2);

                    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot3([pS(1) pD(1)], [pS(2) pD(2)], [pS(3) pD(3)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end
            end

            axis equal
            grid on
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Probe geometry 3D');
            view(3)

        else
            figure('Name','Probe geometry 2D','Color','w'); hold on

            if ~isempty(landPos) && size(landPos,2) >= 2
                scatter(landPos(:,1), landPos(:,2), 30, [0.7 0.7 0.7], 'filled');
            end

            scatter(srcPos(:,1), srcPos(:,2), 90, 'r', 'filled');
            scatter(detPos(:,1), detPos(:,2), 90, 'b', 'filled');

            for i = 1:size(srcPos,1)
                text(srcPos(i,1), srcPos(i,2), sprintf(' S%d', i), ...
                    'Color', 'r', 'FontSize', 9);
            end

            for i = 1:size(detPos,1)
                text(detPos(i,1), detPos(i,2), sprintf(' D%d', i), ...
                    'Color', 'b', 'FontSize', 9);
            end

            if showChannelLinks && ~isempty(ml)
                for ch = 1:size(ml,1)
                    sIdx = ml(ch,1);
                    dIdx = ml(ch,2);

                    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot([pS(1) pD(1)], [pS(2) pD(2)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end
            end

            axis equal
            grid on
            xlabel('X'); ylabel('Y');
            title('Probe geometry 2D');
        end
    end
end

%% ================= PROBE HRF MAP 2D/3D =================
if showProbe3D_HRF
    srcPos  = get_probe_srcPos(probe);
    detPos  = get_probe_detPos(probe);
    landPos = get_probe_landmarkPos(probe);

    yBA_probe  = HRF_by_condition{condIdxPlot};
    tHRF_probe = tHRF_by_condition{condIdxPlot};
    mlBA_probe = mlBA_by_condition{condIdxPlot};

    if isempty(srcPos) || isempty(detPos)
        warning('Probe HRF map skipped: source/detector coordinates not available.');
    elseif isempty(yBA_probe)
        warning('Probe HRF map skipped: HRF not available for selected condition.');
    else
        hbName = hb_idx_to_name(probeMetricHbTypeIdx);

        tMask = tHRF_probe >= probeMetricTimeWindow(1) & tHRF_probe <= probeMetricTimeWindow(2);
        if ~any(tMask)
            warning('No HRF samples in selected time window.');
        else
            nPairs = estimate_num_sd_pairs(mlBA_probe);
            metricAllChannels = nan(1,nPairs);

            for ch = 1:nPairs
                idxHb = find_hb_channels(mlBA_probe, ch, hbName);
                if ~isempty(idxHb)
                    metricAllChannels(ch) = mean(yBA_probe(tMask, idxHb(1)), 1, 'omitnan');
                end
            end

            idxWl = get_meas_idx_for_wavelength(ml, probeUseWavelengthIdx);
            pairMap = unique(ml(:,1:2), 'rows', 'stable');

            % Use all available S-D pairs from BlockAvg output.
            idxUsePairs = 1:nPairs;
            metricUse = metricAllChannels(idxUsePairs);

            [midPos, srcDetPairs] = compute_pair_midpoints(pairMap, srcPos, detPos, idxUsePairs);

            is3D = (size(srcPos,2) >= 3) && (size(detPos,2) >= 3);

            if is3D
                figure('Name','Probe HRF map 3D','Color','w'); hold on

                if ~isempty(landPos) && size(landPos,2) >= 3
                    scatter3(landPos(:,1), landPos(:,2), landPos(:,3), 30, [0.7 0.7 0.7], 'filled');
                end

                scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 70, 'r', 'filled');
                scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 70, 'b', 'filled');

                if showChannelLinks
                    for i = 1:size(srcDetPairs,1)
                        sIdx = srcDetPairs(i,1);
                        dIdx = srcDetPairs(i,2);
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot3([pS(1) pD(1)], [pS(2) pD(2)], [pS(3) pD(3)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end

                scatter3(midPos(:,1), midPos(:,2), midPos(:,3), 120, metricUse(:), 'filled');

                axis equal
                grid on
                xlabel('X'); ylabel('Y'); zlabel('Z');
                title(sprintf('Probe HRF map 3D | %s | %s | mean %.1f-%.1f s', ...
                    condNames(condIdxPlot), hbName, ...
                    probeMetricTimeWindow(1), probeMetricTimeWindow(2)), ...
                    'Interpreter', 'none');
                colorbar
                view(3)

            else
                figure('Name','Probe HRF map 2D','Color','w'); hold on

                if ~isempty(landPos) && size(landPos,2) >= 2
                    scatter(landPos(:,1), landPos(:,2), 30, [0.7 0.7 0.7], 'filled');
                end

                scatter(srcPos(:,1), srcPos(:,2), 70, 'r', 'filled');
                scatter(detPos(:,1), detPos(:,2), 70, 'b', 'filled');

                if showChannelLinks
                    for i = 1:size(srcDetPairs,1)
                        sIdx = srcDetPairs(i,1);
                        dIdx = srcDetPairs(i,2);
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot([pS(1) pD(1)], [pS(2) pD(2)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end

                scatter(midPos(:,1), midPos(:,2), 120, metricUse(:), 'filled');

                axis equal
                grid on
                xlabel('X'); ylabel('Y');
                title(sprintf('Probe HRF map 2D | %s | %s | mean %.1f-%.1f s', ...
                    condNames(condIdxPlot), hbName, ...
                    probeMetricTimeWindow(1), probeMetricTimeWindow(2)), ...
                    'Interpreter', 'none');
                colorbar
            end
        end
    end
end

%% ================= SAVE ALL FIGURES =================
outFigDir = fullfile(pwd, sprintf('Figures_%s', motionCorrectionMethod));
if ~exist(outFigDir, 'dir')
    mkdir(outFigDir);
end

figHandles = findall(0, 'Type', 'figure');

for i = 1:length(figHandles)
    fig = figHandles(i);

    figName = get(fig, 'Name');
    if isempty(figName)
        figName = sprintf('Figure_%d', i);
    end

    % Clean filename for Windows
    figName = char(figName);
    figName = regexprep(figName, '[<>:"/\\|?*]', '_');
    figName = regexprep(figName, '\s+', '_');
    figName = regexprep(figName, '_+', '_');

    % Add index to avoid duplicated names
    figName = sprintf('%02d_%s', i, figName);

    pngFile = fullfile(outFigDir, [figName '.png']);
    figFile = fullfile(outFigDir, [figName '.fig']);

    try
        exportgraphics(fig, pngFile, 'Resolution', 300);
    catch
        saveas(fig, pngFile);
    end

    savefig(fig, figFile);
end

fprintf('All figures saved in folder: %s\n', outFigDir);
%% ================= SAVE =================
save(outMat, ...
    'snirfFile','keepStims', ...
    'useMotionArtifactDetect','useMotionArtifactByChannel', ...
    'motionCorrectionMethod','useBandpassFilter', ...
    'tMotion','tMask','STDEVthresh','AMPthresh', ...
    'iqr_wavelet','turnon_wavelet','pSpline','hpf','lpf','dpf', ...
    'exampleCh','hbTypeIdx','condIdxPlot', ...
    'badChannels','badPairIdx','meanRaw','cvRaw', ...
    'dodRaw','dodWave','dodSpline','dodCorr','dodFilt','ml', ...
    'HRF_by_condition','tHRF_by_condition','mlBA_by_condition', ...
    'Results','SNR_table','probe');

fprintf('Saved: %s\n', outMat);

%% ===================== LOCAL HELPERS =====================

function d0 = get_first_data_element(dataObj)
% Return the first data element robustly.
try
    d0 = dataObj(1);
    return
catch
end

try
    if iscell(dataObj)
        d0 = dataObj{1};
        return
    end
catch
end

d0 = dataObj;
end

function Y = get_data_matrix(dataObj)
% Extract numeric data matrix [time x measurements].
try
    Y = dataObj.GetDataTimeSeries();
    return
catch
end

try
    Y = dataObj.dataTimeSeries;
    return
catch
end

try
    Y = dataObj.GetDataTimeSeries('reshape');
    return
catch
end

error('Unable to extract data matrix from dataObj.');
end

function t = get_time_from_data(dataObj)
% Extract time vector robustly.
try
    t = dataObj.GetTime();
    t = t(:);
    return
catch
end

try
    t = dataObj.time;
    t = t(:);
    return
catch
end

try
    t = dataObj.t;
    t = t(:);
    return
catch
end

error('Unable to extract time vector from dataObj.');
end

function ml = get_meas_list_safe(dataObj)
% Extract measurement list robustly.
try
    ml = dataObj.GetMeasList('reshape');
    return
catch
end

try
    ml = dataObj.GetMeasList();
    return
catch
end

try
    ml = dataObj.measList;
    return
catch
end

error('Unable to extract measurement list from dataObj.');
end

function wl = get_probe_wavelengths(probe)
% Extract wavelength vector from probe.
try
    wl = probe.GetWavelengths();
catch
    try
        wl = probe.wavelengths;
    catch
        try
            wl = probe.lambda;
        catch
            wl = [];
        end
    end
end
wl = wl(:)';
end

function names = get_stim_names(stimArr)
% Return stimulus names as a string array.
if isempty(stimArr)
    names = string.empty(0,1);
    return
end

names = strings(numel(stimArr),1);
for i = 1:numel(stimArr)
    names(i) = get_stim_name_one(stimArr(i));
end
end

function name = get_stim_name_one(stim)
% Extract one stimulus name robustly.
try
    name = string(stim.name);
    return
catch
end

try
    name = string(stim.GetName());
    return
catch
end

name = "stim";
end

function data = get_stim_data_one(stim)
% Extract stimulus data, usually [onset duration amplitude].
try
    data = stim.data;
    return
catch
end

try
    data = stim.GetData();
    return
catch
end

data = [];
end

function ttl = make_channel_title(ml, ch)
% Build a readable title for one measurement.
if isempty(ml)
    ttl = sprintf('Ch %d', ch);
    return
end

if size(ml,2) >= 4
    ttl = sprintf('Ch %d | S%d-D%d | wlIdx %d', ch, ml(ch,1), ml(ch,2), ml(ch,4));
elseif size(ml,2) >= 2
    ttl = sprintf('Ch %d | S%d-D%d', ch, ml(ch,1), ml(ch,2));
else
    ttl = sprintf('Ch %d', ch);
end
end

function srcPos = get_probe_srcPos(probe)
% Extract source positions, preferring 3D then 2D.
srcPos = [];

try
    srcPos = probe.sourcePos3D;
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.GetSrcPos();
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.sourcePos2D;
    if ~isempty(srcPos), return; end
catch
end
end

function detPos = get_probe_detPos(probe)
% Extract detector positions, preferring 3D then 2D.
detPos = [];

try
    detPos = probe.detectorPos3D;
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.GetDetPos();
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.detectorPos2D;
    if ~isempty(detPos), return; end
catch
end
end

function landPos = get_probe_landmarkPos(probe)
% Extract landmark positions, preferring 3D then 2D.
landPos = [];

try
    landPos = probe.landmarkPos3D;
    if ~isempty(landPos), return; end
catch
end

try
    landPos = probe.landmarkPos2D;
    if ~isempty(landPos), return; end
catch
end
end

function idx = get_meas_idx_for_wavelength(ml, wlIdxWanted)
% Return indices of raw OD measurements belonging to a given wavelength index.
idx = [];
if isempty(ml)
    return
end

if size(ml,2) >= 4
    idx = find(ml(:,4) == wlIdxWanted);
end
end

function idx = find_hb_channels(ml, chPairIdx, hbName)
% Find columns corresponding to one S-D pair and one Hb type.
%
% This helper assumes that BlockAvg concentration output is organized so that
% each unique source-detector pair appears up to 3 times in order:
% HbO, HbR, HbT.
%
% This is a practical reconstruction step for Homer3 native BlockAvg output.

idx = [];
if isempty(ml)
    return
end

pairs = unique(ml(:,1:2), 'rows', 'stable');

if chPairIdx > size(pairs,1)
    return
end

targetPair = pairs(chPairIdx,:);
pairIdx = find(ml(:,1)==targetPair(1) & ml(:,2)==targetPair(2));

switch hbName
    case 'HbO'
        if numel(pairIdx) >= 1, idx = pairIdx(1); end
    case 'HbR'
        if numel(pairIdx) >= 2, idx = pairIdx(2); end
    case 'HbT'
        if numel(pairIdx) >= 3, idx = pairIdx(3); end
end
end

function nPairs = estimate_num_sd_pairs(ml)
% Estimate the number of unique source-detector pairs.
if isempty(ml)
    nPairs = 0;
    return
end
pairs = unique(ml(:,1:2), 'rows', 'stable');
nPairs = size(pairs,1);
end

function hbName = hb_idx_to_name(hbTypeIdx)
% Map Hb index to Hb name.
switch hbTypeIdx
    case 1
        hbName = 'HbO';
    case 2
        hbName = 'HbR';
    case 3
        hbName = 'HbT';
    otherwise
        hbName = 'HbO';
end
end

function tIncCh_out = build_spline_tIncCh(dataObj, tInc, tIncCh)
% Build a channel-wise inclusion mask in the format expected by Homer3.
if ~isempty(tIncCh)
    if iscell(tIncCh)
        tIncCh_out = tIncCh;
    else
        tIncCh_out = {tIncCh};
    end
    return
end

if isempty(tInc)
    tIncCh_out = [];
    return
end

if ~iscell(tInc)
    tInc = {tInc};
end

tIncCh_out = cell(length(dataObj), 1);

for iBlk = 1:length(dataObj)
    d = dataObj(iBlk).GetDataTimeSeries();
    nTime = size(d, 1);
    nMeas = size(d, 2);

    if iBlk <= numel(tInc) && ~isempty(tInc{iBlk})
        tMask = tInc{iBlk};
    else
        tMask = ones(nTime, 1);
    end

    tMask = tMask(:);
    if numel(tMask) < nTime
        tMask(end+1:nTime, 1) = 1;
    elseif numel(tMask) > nTime
        tMask = tMask(1:nTime);
    end

    tIncCh_out{iBlk} = repmat(tMask, 1, nMeas);
end
end

function [midPos, srcDetPairs] = compute_pair_midpoints(pairMap, srcPos, detPos, idxPairs)
% Compute midpoints for source-detector pairs.
n = numel(idxPairs);
midPos = nan(n,3);
srcDetPairs = nan(n,2);

for i = 1:n
    pairIdx = idxPairs(i);

    if pairIdx > size(pairMap,1)
        continue
    end

    sIdx = pairMap(pairIdx,1);
    dIdx = pairMap(pairIdx,2);

    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
        pS = srcPos(sIdx,:);
        pD = detPos(dIdx,:);

        if size(pS,2) == 2
            pS(3) = 0;
        end
        if size(pD,2) == 2
            pD(3) = 0;
        end

        midPos(i,:) = (pS + pD)/2;
        srcDetPairs(i,:) = [sIdx dIdx];
    end
end
end
