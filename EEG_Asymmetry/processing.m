% =========================================================================
% EEG Batch Preprocessing Script (Group D) - Compatible BIDS File Loop
% =========================================================================

clear; clc;

% Add vendored EEGLAB to path if needed
if isempty(which('eeglab'))
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'eeglab2026.0.0'));
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % Initialize EEGLAB

% Define directories (relative to this script location)
base_dir = fileparts(mfilename('fullpath'));
raw_data_dir = fullfile(base_dir, 'DATA');
processed_dir = fullfile(base_dir, 'PROCESSED_DATA');

if ~exist(processed_dir, 'dir')
    mkdir(processed_dir);
end

% Recursively search through all subfolders for .vhdr files
all_files = dir(fullfile(raw_data_dir, '**', '*.vhdr'));
num_files = length(all_files);

% Start the For Loop
for f = 1:num_files
    
    current_file = all_files(f).name;
    current_folder = all_files(f).folder;
    
    % --- FILTERING LOGIC ---
    % 1. Check if the file is one of your target tasks
    is_target_task = contains(current_file, 'task-eyesclosed') || contains(current_file, 'task-mathematic');
    
    % 2. Check if the folder is from session 1
    is_target_session = contains(current_folder, 'ses-session1');
    
    % If it is NOT the right task OR NOT the right session, skip it
    if ~(is_target_task && is_target_session)
        continue; 
    end
    % -----------------------
    
    % Extract metadata parts from the BIDS filename for clean reporting
    name_parts = strsplit(current_file, '_');
    subject_id = name_parts{1}; % 'sub-01'
    session_id = name_parts{2}; % 'ses-session1'
    task_id    = name_parts{3}; % 'task-eyesclosed' or 'task-mathematic'
    
    % Create clean output filename
    output_filename = sprintf('%s_%s_%s_processed.set', subject_id, session_id, task_id);
    
    fprintf('\n==================================================\n');
    fprintf('Processing File %d of %d\n', f, num_files);
    fprintf('Subject: %s | Session: %s | %s\n', subject_id, session_id, task_id);
    fprintf('==================================================\n');
    
    % --- STEP 1: Import Data using BrainVision native loader ---
    EEG = pop_loadbv(current_folder, current_file);
    
    % --- STEP 2: Dynamically load the associated BIDS electrodes .tsv ---
    tsv_filename = sprintf('%s_%s_electrodes.tsv', subject_id, session_id);
    tsv_filepath = fullfile(current_folder, tsv_filename);
    
    if isfile(tsv_filepath)
        fprintf('  -> Applying exact BIDS channel locations from: %s\n', tsv_filename);
        EEG = pop_chanedit(EEG, 'load', {tsv_filepath, 'filetype', 'autodetect'});
    else
        fprintf('  -> WARNING: %s not found. Falling back to generic template.\n', tsv_filename);
        EEG = pop_chanedit(EEG, 'lookup','standard_1005.elc');
    end
    
    % --- STEP 3 & 4: Filtering ---
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 70);
    
    % --- STEP 5: Resampling ---
    EEG = pop_resample(EEG, 250);
    
    % --- STEP 6: Artifact Removal using ASR ---
    orig_chanlocs = EEG.chanlocs;
    orig_chans = {EEG.chanlocs.labels};
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, ...
        'LineNoiseCriterion', 4, 'Highpass', 'off', 'BurstCriterion', 20, ...
        'WindowCriterion', 'off', 'BurstRejection', 'on', 'Distance', 'Euclidian');
    
    clean_chans = {EEG.chanlocs.labels};
    removed_chans = setdiff(orig_chans, clean_chans);
    
    if isempty(removed_chans)
        fprintf('  -> ASR: No bad channels were removed.\n');
    else
        fprintf('  -> ASR WARNING: Removed %d bad channel(s): %s\n', ...
            length(removed_chans), strjoin(removed_chans, ', '));
    end
    EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
    fprintf('  -> Interpolated channels to full montage (%d channels).\n', length(orig_chanlocs));
    
    % --- STEP 7: ICA ---
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
    
    % --- STEP 8: ICLabel & Component Rejection ---
    EEG = pop_iclabel(EEG, 'default');
    thresholds = [NaN NaN; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; NaN NaN];
    EEG = pop_icflag(EEG, thresholds);
    rej_comps = find(EEG.reject.gcompreject);
    
    if isempty(rej_comps)
        fprintf('  -> ICA: No noisy components identified.\n');
    else
        fprintf('  -> ICA: Removed %d noise component(s): %s\n', ...
            length(rej_comps), num2str(rej_comps));
    end
    EEG = pop_subcomp(EEG, rej_comps, 0);
    
    % --- STEP 9: Re-referencing ---
    EEG = pop_reref(EEG, []);
    
    % --- STEP 10: Save Processed Dataset ---
    EEG = pop_saveset(EEG, 'filename', output_filename, 'filepath', processed_dir);
    
    % EEG is overwritten each iteration; no need to manually clear ALLEEG
end

fprintf('\nBIDS-compatible batch preprocessing completed successfully!\n');
