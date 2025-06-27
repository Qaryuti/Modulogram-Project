% batch_modulogram.m

% Define parameters
subjects = {'EMU024', 'EMU001', 'EMU025', 'EMU030'};
alignments = {'win', 'loss'};
modulatedRanges = {[30 120], [20 80]};
sessionNums = 1:3;
bandwidth = 5;
numBins = 36;

% Base directories
baseResultsDir = '/home/alabwam1/Desktop/modulogram_V7/trial-run-3';
dataDir = '/media/Data/Human_Intracranial_MAD';

% Optional: open a master log file
masterLogPath = fullfile(baseResultsDir, 'batch_log.txt');
masterLog = fopen(masterLogPath, 'a');
fprintf(masterLog, '\n==== Batch run started: %s ====\n', datestr(now));

% Loop over all combinations
for s = 1:length(subjects)
    for a = 1:length(alignments)
        for m = 1:length(modulatedRanges)
            for sess = sessionNums

                subject = subjects{s};
                alignment = alignments{a};
                modRange = modulatedRanges{m};

                % Descriptive results folder name
                resultFolderName = sprintf('%s_%s_mod%d-%d_session%d_bw%d_bins%d', ...
                    subject, alignment, modRange(1), modRange(2), sess, bandwidth, numBins);
                resultsDir = fullfile(baseResultsDir, resultFolderName);

                % Create output directory if it doesn't exist
                if ~exist(resultsDir, 'dir')
                    mkdir(resultsDir);
                end

                % Config struct
                config.subjects = {subject};
                config.sessionNum = sess;
                config.alignments = {alignment};
                config.dataDir = dataDir;
                config.reference = 'Ground';
                config.resultsDir = resultsDir;
                config.modulatedRange = modRange;
                config.modWaveRanges.Theta = [4 8];
                config.bandwidth = bandwidth;
                config.numBins = numBins;
                config.windowParams = struct('window', [0 1], 'excludeOutOfBounds', true);

                % Log path for this specific run
                localLogPath = fullfile(resultsDir, 'run_log.txt');
                localLog = fopen(localLogPath, 'a');

                % Log configuration
                configText = sprintf([ ...
                    '\n--- Run started at %s ---\n', ...
                    'Subject: %s\nAlignment: %s\nSession: %d\nModulated Range: [%d %d]\n', ...
                    'Bandwidth: %d\nBins: %d\nOutput Dir: %s\n'], ...
                    datestr(now), subject, alignment, sess, ...
                    modRange(1), modRange(2), bandwidth, numBins, resultsDir);

                fprintf(configText);
                fprintf(localLog, configText);
                fprintf(masterLog, configText);

                % Try-catch to handle errors
                try
                    modulogram_pipeline(config);
                    fprintf(localLog, '✓ Completed successfully at %s\n', datestr(now));
                    fprintf(masterLog, '✓ Completed: %s at %s\n', resultFolderName, datestr(now));
                catch ME
                    fprintf(localLog, '✗ ERROR at %s:\n%s\n', datestr(now), getReport(ME));
                    fprintf(masterLog, '✗ ERROR in %s at %s\n%s\n', resultFolderName, datestr(now), getReport(ME));
                end

                fclose(localLog);
            end
        end
    end
end

% Close the master log
fprintf(masterLog, '==== Batch run ended: %s ====\n\n', datestr(now));
fclose(masterLog);
