config = struct();
config.subjectList = {'EMU024', 'EMU025', 'EMU030'};
config.sessionNums = [1];  % Or [1 2 3] if needed
config.lossRoot    = 'loss-1';         % This folder is inside modulogram_V6
config.winRoot     = 'win-1';
config.outputRoot  = 'results_stats';  % Create a new folder if you want


run_pac_stats_all(config);

function run_pac_stats_all(config)
    % Loop over all subjects and sessions
    for subjIdx = 1:numel(config.subjectList)
        subjectID = config.subjectList{subjIdx};

        for s = 1:numel(config.sessionNums)
            sesnum = config.sessionNums(s);

            % Build loss/win paths
            lossFile = fullfile(config.lossRoot, sprintf('%s_session%d_struct.mat', subjectID, sesnum));
            winFile  = fullfile(config.winRoot,  sprintf('%s_session%d_struct.mat', subjectID, sesnum));

            % Check existence
            if ~exist(lossFile, 'file')
                fprintf('Missing loss file for %s, session %d\n', subjectID, sesnum);
                continue;
            end
            if ~exist(winFile, 'file')
                fprintf('Missing win file for %s, session %d\n', subjectID, sesnum);
                continue;
            end

            % Create inner config and run
            innerConfig = struct();
            innerConfig.lossPath    = lossFile;
            innerConfig.winPath     = winFile;
            innerConfig.subjectID   = subjectID;
            innerConfig.sessionNum  = sesnum;
            innerConfig.outputRoot  = config.outputRoot;

            fprintf('Running PAC stats for %s session %d...\n', subjectID, sesnum);
            run_pac_stats(innerConfig);
        end
    end
end
