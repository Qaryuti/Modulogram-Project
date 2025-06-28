function outFile = convert_modulogram_for_logistic(modulogramMatFile)
% CONVERT_MODULOGRAM_FOR_LOGISTIC
%   Converts a nested modulogram .mat file into the simpler
%   sessStruct format expected by analyze_winloss_logistic_pipeline,
%   using folder name to assign labels.
%
%   outFile = convert_modulogram_for_logistic(modulogramMatFile)
%
%   Saves the result under ../Data relative to this script.
%
%   Example:
%       outFile = convert_modulogram_for_logistic('/path/to/EMU024_session1_struct.mat');

% Load the nested modulogram structure
fprintf('Loading modulogram data from %s\n', modulogramMatFile);
data = load(modulogramMatFile);

if ~isfield(data,'allData')
    error('File does not contain allData');
end
allData = data.allData;

% Parse subject/session/alignment from the filename
[~,baseName] = fileparts(modulogramMatFile);
tokens = regexp(baseName, '(EMU\d+)_session(\d+)', 'tokens', 'once');
if isempty(tokens)
    error('Could not parse subjectID/sessionNum from filename');
end
subjectID = tokens{1};
sessionNum = str2double(tokens{2});

% determine alignment type from the path
if contains(modulogramMatFile,'win','IgnoreCase',true)
    alignment = 'win';
    labelValue = 1;
elseif contains(modulogramMatFile,'loss','IgnoreCase',true)
    alignment = 'loss';
    labelValue = 0;
else
    error('Cannot infer alignment (win/loss) from file path');
end

% initialize
sessStruct.subjectID = subjectID;
sessStruct.sessionNum = sessionNum;
sessStruct.channels = {};

% pull channels
chanStructs = allData.(subjectID).session(sessionNum).alignment.(alignment).channel;
chanNames = fieldnames(chanStructs);

for c = 1:numel(chanNames)
    chName = chanNames{c};
    ch = chanStructs.(chName);
    
    channel.finalAggMI = ch.finalAggMI;
    channel.finalAggAmplP = ch.finalAggAmplP;
    channel.anatomicalRegion = ch.anatomicalRegion;
    channel.trialMI = ch.MI_per_trial';
    
    nTrials = size(channel.trialMI, 1);
    channel.trialLabels = repmat(labelValue, nTrials, 1);
    
    sessStruct.channels{end+1} = channel;
end

% Save into Data folder
analysisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(analysisDir,'Data');
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end
outFile = fullfile(dataDir, sprintf('%s_session%d_sessStruct.mat',subjectID,sessionNum));
save(outFile,'sessStruct');
fprintf('Saved converted sessStruct to %s\n', outFile);

end
