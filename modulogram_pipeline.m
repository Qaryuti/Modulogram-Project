function modulogram_pipeline(config)
%MODULOGRAM_PIPELINE Orchestrate modulogram computation for one session.
%   MODULOGRAM_PIPELINE(CONFIG) loops over subjects and channels for a
%   specified session and alignments. Helper functions live in src.

% Ensure helper functions in the src directory are on the MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'src'));

% Default FIR filter order if not specified
if ~isfield(config, 'fir_order')
    config.fir_order = 1000;
end

% Enable optional debug mode for verbose logging and plotting
if ~isfield(config, 'debug_mode')
    config.debug_mode = false;
end

if ~isfield(config, 'logFile')
    ts = datestr(now,'yyyymmdd_HHMMSS');
    if ~exist(config.resultsDir,'dir'); mkdir(config.resultsDir); end
    config.logFile = fullfile(config.resultsDir, ['pipeline_' ts '.log']);
end
logFID = fopen(config.logFile,'w');

logf('===== Modulogram Pipeline =====\n');
logf('Subjects: %d\n', numel(config.subjects));
logf('Session : %d\n', config.sessionNum);
logf('Alignments: %s\n', strjoin(config.alignments, ', '));
if config.debug_mode
    logf('Debug mode ENABLED - real time figures will be displayed.\n');
end

subjectParams = struct( ...
    'EMU001', struct('num_sessions',3,'Fs',1000), ...
    'EMU024', struct('num_sessions',3,'Fs',2048), ...
    'EMU025', struct('num_sessions',2,'Fs',2048), ...
    'EMU030', struct('num_sessions',2,'Fs',2048), ...
    'EMU036', struct('num_sessions',4,'Fs',2048), ...
    'EMU037', struct('num_sessions',4,'Fs',2048), ...
    'EMU038', struct('num_sessions',1,'Fs',2048), ...
    'EMU039', struct('num_sessions',4,'Fs',2048), ...
    'EMU040', struct('num_sessions',2,'Fs',2048), ...
    'EMU041', struct('num_sessions',9,'Fs',2048), ...
    'EMU047', struct('num_sessions',1,'Fs',2048), ...
    'EMU051', struct('num_sessions',1,'Fs',2048) ...
    );

for subjIdx = 1:numel(config.subjects)
    subjectID = config.subjects{subjIdx};
    if ~isfield(subjectParams, subjectID)
        logf('No parameters defined for subject %s, skipping.\n', subjectID);
        continue;
    end
    Fs = subjectParams.(subjectID).Fs;
    sesnum = config.sessionNum;
    logf('\n=== [%d/%d] Subject %s | Session %d ===\n', subjIdx, numel(config.subjects), subjectID, sesnum);

    data_base_dir = fullfile(config.dataDir, '1_formatted');
    if strcmpi(config.reference, 'Ground')
        setupFile = fullfile(data_base_dir, subjectID, sprintf('%s_MAD_SES%d_Setup.mat', subjectID, sesnum));
    else
        setupFile = fullfile(data_base_dir, subjectID, sprintf('%s_MAD_SES%d_Setup_%s.mat', subjectID, sesnum, config.reference));
    end
    if ~exist(setupFile, 'file')
        logf('    Setup file missing for %s, Session %d\n', subjectID, sesnum);
        continue;
    end
    load(setupFile, 'filters', 'trial_times', 'trial_words', 'elec_area', 'elec_ind');
    logf('  Loaded setup: %d channels\n', numel(elec_ind));
    allData.(subjectID).metadata = struct('Fs', Fs, 'subjectID', subjectID);
    allData.(subjectID).session(sesnum).metadata = struct( ...
        'elec_area', elec_area, 'elec_ind', elec_ind, ...
        'filters', filters, 'trial_times', {trial_times}, 'trial_words', {trial_words});

    for alignIdx = 1:numel(config.alignments)
        alignment = config.alignments{alignIdx};
        logf('  -> Alignment %d/%d: %s\n', alignIdx, numel(config.alignments), alignment);
        for chIdx = 1:numel(elec_ind)
            channelLabel = sprintf('%03d', elec_ind(chIdx));
            if exist('elec_area','var') && ~isempty(elec_area)
                if iscell(elec_area)
                    anatomicalRegion = elec_area{chIdx};
                else
                    anatomicalRegion = strtrim(elec_area(chIdx,:));
                end
            else
                anatomicalRegion = 'Unknown';
            end

            logf('    -> Channel %s (%s) [%d/%d]\n', channelLabel, anatomicalRegion, chIdx, numel(elec_ind));
            modStruct = create_single_modulogram(config, subjectID, sesnum, channelLabel, alignment, Fs, filters, trial_times, trial_words, anatomicalRegion);
            if isempty(modStruct)
                logf('       [!] No valid trials for %s, Session %d, Channel %s\n', subjectID, sesnum, channelLabel);
                continue;
            end
            chanField = ['CH' channelLabel];
            allData.(subjectID).session(sesnum).alignment.(alignment).channel.(chanField) = modStruct;
        end
    end
    subjectFile = fullfile(config.resultsDir, sprintf('%s_session%d_struct.mat', subjectID, sesnum));
    save(subjectFile, 'allData');
    logf('  Saved subject data to %s\n', subjectFile);
end
tsFinal = datestr(now,'yyyymmdd_HHMMSS');
finalFile = fullfile(config.resultsDir, sprintf('allSubjects_session%d_%s.mat', config.sessionNum, tsFinal));
save(finalFile, 'allData');
logf('Saved combined data to %s\n', finalFile);
fclose(logFID);

    function logf(fmt, varargin)
        timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
        msg = sprintf(fmt, varargin{:});
        fprintf('%s', ['[' timestamp '] ' msg]);
        if ~isempty(logFID) && logFID>0
            fprintf(logFID, '%s', ['[' timestamp '] ' msg]);
        end
    end
end
