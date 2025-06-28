function outFile = convert_modulogram_merge(winMatFile, lossMatFile)
% CONVERT_MODULOGRAM_MERGE
%   Combines win-aligned and loss-aligned modulogram results
%   into a single sessStruct for logistic regression, labeling
%   purely based on folder name.
%
%   outFile = convert_modulogram_merge(winMatFile, lossMatFile)
%
%   Saves under analysis/Data

fprintf('\n=== Starting MERGE ===\n');
fprintf('Loading WIN modulogram: %s\n', winMatFile);
dataWin = load(winMatFile);
fprintf('Loading LOSS modulogram: %s\n', lossMatFile);
dataLoss = load(lossMatFile);

if ~isfield(dataWin,'allData') || ~isfield(dataLoss,'allData')
    error('Missing allData in one of the inputs.');
end

% parse subject/session
[~,baseName] = fileparts(winMatFile);
tokens = regexp(baseName, '(EMU\d+)_session(\d+)', 'tokens', 'once');
if isempty(tokens)
    error('Could not parse subjectID/sessionNum');
end
subjectID = tokens{1};
sessionNum = str2double(tokens{2});

sessStruct.subjectID = subjectID;
sessStruct.sessionNum = sessionNum;
sessStruct.channels = {};

% get channels
winChannels = dataWin.allData.(subjectID).session(sessionNum).alignment.win.channel;
lossChannels = dataLoss.allData.(subjectID).session(sessionNum).alignment.loss.channel;

allChanNames = union(fieldnames(winChannels), fieldnames(lossChannels));
fprintf('Merging %d unique channels for session %d\n', numel(allChanNames), sessionNum);

for c = 1:numel(allChanNames)
    chName = allChanNames{c};
    fprintf('\n-- Channel: %s --\n', chName);
    
    channel.trialMI = [];
    channel.trialLabels = [];
    channel.finalAggMI = [];
    channel.finalAggAmplP = [];
    channel.anatomicalRegion = '';

    % merge WIN trials
    if isfield(winChannels, chName)
        ch = winChannels.(chName);
        mi = ch.MI_per_trial';
        nTrials = size(mi,1);
        labels = ones(nTrials,1); % assign 1
        
        fprintf('  WIN trials found: %d (assigned label 1)\n', nTrials);
        fprintf('    trialNums: %s\n', mat2str(ch.trialNumbers));
        
        channel.trialMI = [channel.trialMI; mi];
        channel.trialLabels = [channel.trialLabels; labels];
        
        channel.finalAggMI = ch.finalAggMI;
        channel.finalAggAmplP = ch.finalAggAmplP;
        channel.anatomicalRegion = ch.anatomicalRegion;
    else
        fprintf('  WIN trials missing for channel %s\n', chName);
    end
    
    % merge LOSS trials
    if isfield(lossChannels, chName)
        ch = lossChannels.(chName);
        mi = ch.MI_per_trial';
        nTrials = size(mi,1);
        labels = zeros(nTrials,1); % assign 0
        
        fprintf('  LOSS trials found: %d (assigned label 0)\n', nTrials);
        fprintf('    trialNums: %s\n', mat2str(ch.trialNumbers));
        
        channel.trialMI = [channel.trialMI; mi];
        channel.trialLabels = [channel.trialLabels; labels];
        
        % preserve metadata if needed
        if isempty(channel.finalAggMI)
            channel.finalAggMI = ch.finalAggMI;
            channel.finalAggAmplP = ch.finalAggAmplP;
            channel.anatomicalRegion = ch.anatomicalRegion;
        end
    else
        fprintf('  LOSS trials missing for channel %s\n', chName);
    end
    
    % summary
    fprintf('  Combined trialLabels in channel: %s\n', mat2str(channel.trialLabels));
    fprintf('  Unique labels: %s\n', mat2str(unique(channel.trialLabels)));
    
    sessStruct.channels{end+1} = channel;
end

% save
dataDir = fullfile(fileparts(mfilename('fullpath')), 'Data');
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end
outFile = fullfile(dataDir, sprintf('%s_session%d_sessStruct.mat',subjectID,sessionNum));
save(outFile, 'sessStruct');
fprintf('Saved combined sessStruct to %s\n', outFile);
fprintf('=== MERGE COMPLETED ===\n');

end
