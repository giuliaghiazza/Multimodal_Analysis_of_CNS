%% Simple SNIRF viewer with Homer3
% Load one SNIRF file with Homer3 and visualize the raw signals.

clear; close all; clc

%% ================= PATHS =================
HOMER3_DIR = '/home/giulia/MATLAB/Toolboxes/Homer3';
PATCH_DIR = fullfile(fileparts(mfilename('fullpath')), 'matlab_compat_patches');
if exist(PATCH_DIR, 'dir') == 7
    addpath(PATCH_DIR, '-begin');
end
addpath(genpath(HOMER3_DIR));

%% ================= INPUT =================
snirfFile = fullfile(fileparts(mfilename('fullpath')), 'CW - Motor', '2025-04-01_001.snirf');
exampleMeas = 1;

assert(exist(snirfFile, 'file') == 2, 'SNIRF not found: %s', snirfFile);

%% ================= LOAD SNIRF =================
snirf = SnirfClass(snirfFile);
dataObj = snirf.data;
probe = snirf.probe;

d0 = get_first_data_element(dataObj);
t = get_time_from_data(d0);
Y = get_data_matrix(d0);
ml = get_meas_list_safe(d0);

nTime = size(Y, 1);
nMeas = size(Y, 2);

assert(exampleMeas >= 1 && exampleMeas <= nMeas, ...
    'exampleMeas must be between 1 and %d.', nMeas);

fprintf('Loaded file: %s\n', snirfFile);
fprintf('Time points: %d\n', nTime);
fprintf('Measurements: %d\n', nMeas);

try
    stimNames = strings(numel(snirf.stim), 1);
    for i = 1:numel(snirf.stim)
        stimNames(i) = string(snirf.stim(i).GetName());
    end
    fprintf('Stim names:\n');
    disp(stimNames)
catch
    fprintf('No readable stimulus information found.\n');
end

%% ================= PLOTS =================
figure('Name', 'SNIRF raw signals - all measurements', 'Color', 'w');
plot(t, Y, 'Color', [0.75 0.75 0.75]);
grid on
xlabel('Time (s)');
ylabel('Intensity');
title('Raw SNIRF signals - all measurements');

figure('Name', sprintf('SNIRF raw signal - measurement %d', exampleMeas), 'Color', 'w');
plot(t, Y(:, exampleMeas), 'b-', 'LineWidth', 1.2);
grid on
xlabel('Time (s)');
ylabel('Intensity');
title(sprintf('Raw signal | %s', make_measurement_title(ml, exampleMeas)), ...
    'Interpreter', 'none');

% Simple heatmap for a quick overview of all channels
figure('Name', 'SNIRF raw signals heatmap', 'Color', 'w');
imagesc(t, 1:nMeas, Y');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Raw SNIRF signals heatmap');
colorbar

% Optional probe geometry if positions are available
srcPos = get_probe_srcPos(probe);
detPos = get_probe_detPos(probe);

if ~isempty(srcPos) && ~isempty(detPos)
    figure('Name', 'Probe geometry', 'Color', 'w');
    hold on

    if size(srcPos, 2) >= 3 && size(detPos, 2) >= 3
        scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 90, 'r', 'filled');
        scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 90, 'b', 'filled');

        for i = 1:size(srcPos, 1)
            text(srcPos(i,1), srcPos(i,2), srcPos(i,3), sprintf(' S%d', i), ...
                'Color', 'r');
        end
        for i = 1:size(detPos, 1)
            text(detPos(i,1), detPos(i,2), detPos(i,3), sprintf(' D%d', i), ...
                'Color', 'b');
        end

        for i = 1:size(ml, 1)
            sIdx = ml(i,1);
            dIdx = ml(i,2);
            if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                pS = srcPos(sIdx,:);
                pD = detPos(dIdx,:);
                plot3([pS(1) pD(1)], [pS(2) pD(2)], [pS(3) pD(3)], ...
                    '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
            end
        end

        axis equal
        view(3)
        zlabel('Z');
    else
        scatter(srcPos(:,1), srcPos(:,2), 90, 'r', 'filled');
        scatter(detPos(:,1), detPos(:,2), 90, 'b', 'filled');

        for i = 1:size(srcPos, 1)
            text(srcPos(i,1), srcPos(i,2), sprintf(' S%d', i), 'Color', 'r');
        end
        for i = 1:size(detPos, 1)
            text(detPos(i,1), detPos(i,2), sprintf(' D%d', i), 'Color', 'b');
        end

        for i = 1:size(ml, 1)
            sIdx = ml(i,1);
            dIdx = ml(i,2);
            if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                pS = srcPos(sIdx,:);
                pD = detPos(dIdx,:);
                plot([pS(1) pD(1)], [pS(2) pD(2)], ...
                    '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
            end
        end

        axis equal
    end

    grid on
    xlabel('X');
    ylabel('Y');
    title('Probe geometry');
end

%% ================= LOCAL HELPERS =================
function d0 = get_first_data_element(dataObj)
try
    d0 = dataObj(1);
    return
catch
end

try
    d0 = dataObj.GetData(1);
    return
catch
end

error('Unable to extract the first data object from SNIRF.');
end

function Y = get_data_matrix(dataObj)
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

error('Unable to extract data matrix from SNIRF data object.');
end

function t = get_time_from_data(dataObj)
try
    t = dataObj.GetTime();
    return
catch
end

try
    t = dataObj.time;
    return
catch
end

error('Unable to extract time vector from SNIRF data object.');
end

function ml = get_meas_list_safe(dataObj)
try
    ml = dataObj.GetMeasList();
    return
catch
end

try
    ml = dataObj.measurementList;
    return
catch
end

error('Unable to extract measurement list from SNIRF data object.');
end

function ttl = make_measurement_title(ml, idx)
if isempty(ml)
    ttl = sprintf('Measurement %d', idx);
    return
end

if size(ml, 2) >= 4
    ttl = sprintf('Measurement %d | S%d-D%d | wlIdx %d', ...
        idx, ml(idx,1), ml(idx,2), ml(idx,4));
else
    ttl = sprintf('Measurement %d | S%d-D%d', idx, ml(idx,1), ml(idx,2));
end
end

function srcPos = get_probe_srcPos(probe)
srcPos = [];

try
    srcPos = probe.sourcePos3D;
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.sourcePos2D;
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.GetSrcPos();
catch
end
end

function detPos = get_probe_detPos(probe)
detPos = [];

try
    detPos = probe.detectorPos3D;
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.detectorPos2D;
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.GetDetPos();
catch
end
end
