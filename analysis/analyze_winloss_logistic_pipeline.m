function analyze_winloss_logistic_pipeline(sessionStructPath, resultsDir)
%ANALYZE_WINLOSS_LOGISTIC_PIPELINE Run win/loss logistic regression on MI data.
%   ANALYZE_WINLOSS_LOGISTIC_PIPELINE(SESSIONSTRUCTPATH, RESULTSDIR) loads a
%   session structure containing trial-level modulation index (MI) data and
%   fits logistic regression models predicting win/loss outcomes from MI.of
%
%   Results are saved under RESULTSDIR/Session_X/ with subfolders for all
%   channels ("ALL") and orbitofrontal cortex channels ("OFC"). The function
%   writes model coefficients, ROC curves, confusion matrices, descriptive
%   statistics and a summary log. Extensive debug information is printed
%   during processing.
%
%   The expected input MAT-file should contain a structure with fields:
%       - subjectID      : string identifying the subject
%       - sessionNum     : numeric session number
%       - channels       : cell array of channel-level structs. Each struct
%                          may contain:
%             * finalAggMI       - [1 x numGammaBands] average MI
%             * finalAggAmplP    - [numGammaBands x numPhaseBins]
%             * anatomicalRegion - char or string label
%             * trialMI          - [numTrials x numGammaBands] (optional)
%             * trialLabels      - [numTrials x 1] vector, 1=win, 0=loss
%
%   Example:
%       analyze_winloss_logistic_pipeline('EMU024/Session_1/EMU024_session1_struct.mat', ...
%                                         fullfile('analysis','results'))
%
%   This function requires the Statistics and Machine Learning Toolbox.

if nargin < 2
    error('Two input arguments required: sessionStructPath and resultsDir');
end

%% Load session structure
fprintf('Loading session file: %s\n', sessionStructPath);
loaded = load(sessionStructPath);
flds = fieldnames(loaded);
if numel(flds)==1
    sessStruct = loaded.(flds{1});
else
    sessStruct = loaded;
end

%% Attempt to parse subject and session from path if not present
[~,baseName] = fileparts(sessionStructPath);
tokens = regexp(baseName, '(EMU\d+)_session(\d+)', 'tokens', 'once');
if isempty(tokens)
    tokens = regexp(sessionStructPath, '(EMU\d+).*session(\d+)', 'tokens', 'once');
end
if isempty(tokens)
    error('Unable to parse subjectID and session number from path.');
end
subjectID = tokens{1};
sessionNum = str2double(tokens{2});

% Override if fields exist in structure
if isfield(sessStruct, 'subjectID'); subjectID = sessStruct.subjectID; end
if isfield(sessStruct, 'sessionNum'); sessionNum = sessStruct.sessionNum; end

fprintf('Subject ID : %s\n', subjectID);
fprintf('Session Num: %d\n', sessionNum);

if ~isfield(sessStruct,'channels')
    error('Input structure does not contain a ''channels'' field');
end
channels = sessStruct.channels;
if ~iscell(channels)
    channels = {channels};
end

ofcRegex = '(?i)ofc|orbitofrontal';
decisionRegex = '(?i)ofc|orbitofrontal|insula|amygdala|posterior cingulate|ventral insula';

allIdx = true(1,numel(channels));
ofcIdx = false(1,numel(channels));
decisionIdx = false(1,numel(channels));

for c = 1:numel(channels)
    region = '';
    try
        region = char(channels{c}.anatomicalRegion);
    catch
    end
    if ~isempty(regexp(region, ofcRegex, 'once'))
        ofcIdx(c) = true;
    end
    if ~isempty(regexp(region, decisionRegex, 'once'))
        decisionIdx(c) = true;
    end
end

levels = {
    struct('name','ALL','mask',allIdx), 
    struct('name','OFC','mask',ofcIdx),
    struct('name','DECISION','mask',decisionIdx)
};

for lvl = 1:numel(levels)
    lvlName = levels{lvl}.name;
    chanMask = levels{lvl}.mask;
    fprintf('\n=== Analyzing %s channels (%s) ===\n', lvlName, lvlName);
    chList = channels(chanMask);
    fprintf('Number of channels pooled: %d\n', numel(chList));

    pooledTrialMI = [];
    pooledTrialLabels = [];
    pooledRegions = {};

    for c = 1:numel(chList)
        chan = chList{c};
        try
            mi = chan.trialMI;
            labels = chan.trialLabels;
            pooledTrialMI = [pooledTrialMI; mi];
            pooledTrialLabels = [pooledTrialLabels; labels(:)];
            pooledRegions{end+1} = char(chan.anatomicalRegion); %#ok<AGROW>
        catch ME
            fprintf('Channel %d missing trial data: %s\n', c, ME.message);
        end
    end

    fprintf('Total trials pooled: %d\n', size(pooledTrialMI,1));
    if isempty(pooledTrialMI)
        fprintf('No trial-level data available for %s. Skipping.\n', lvlName);
        continue;
    end

    %% Logistic regression
    nBands = size(pooledTrialMI,2);
    predNames = strcat('Band', string(1:nBands));
    tbl = array2table(pooledTrialMI, 'VariableNames', predNames);
    tbl.Label = pooledTrialLabels;

    formula = 'Label ~ 1';
    for b = 1:nBands
        formula = [formula ' + ' predNames{b}]; %#ok<AGROW>
    end
    fprintf('Fitting logistic regression with predictors: %s\n', strjoin(predNames, ', '));
    mdl = fitglm(tbl, formula, 'Distribution','binomial');

    %% Result directory
    sessDir = fullfile(resultsDir, sprintf('Session_%d', sessionNum), lvlName);
    if ~exist(sessDir,'dir'); mkdir(sessDir); end

    %% Save model
    save(fullfile(sessDir,'logistic_model.mat'), 'mdl');

    %% Save coefficients
    coefTbl = mdl.Coefficients;
    coefFile = fullfile(sessDir, 'coefficients.csv');
    writetable(coefTbl, coefFile);

    %% ROC curve
    [X,Y,~,AUC] = perfcurve(tbl.Label, mdl.Fitted.Probability, 1);
    fig = figure('Visible','off');
    plot(X,Y,'LineWidth',2); grid on;
    xlabel('False positive rate'); ylabel('True positive rate');
    title(sprintf('ROC Curve (AUC = %.3f) - %s', AUC, lvlName));
    rocPath = fullfile(sessDir,'roc_curve.png');
    print(fig, rocPath,'-dpng'); close(fig);
    fprintf('Saved ROC curve to %s\n', rocPath);

    %% Confusion matrix
    predLabels = double(mdl.Fitted.Probability >= 0.5);
    cm = confusionmat(tbl.Label, predLabels);
    fig = figure('Visible','off');
    confusionchart(cm, {'Loss','Win'});
    title(sprintf('Confusion Matrix - %s', lvlName));
    cmPath = fullfile(sessDir,'confusion_matrix.png');
    print(fig, cmPath,'-dpng'); close(fig);
    fprintf('Saved confusion matrix to %s\n', cmPath);

    %% Descriptive statistics
    winIdx = tbl.Label==1;
    lossIdx = tbl.Label==0;
    winMI = pooledTrialMI(winIdx,:);
    lossMI = pooledTrialMI(lossIdx,:);

stats = table((1:nBands)', ...
    mean(pooledTrialMI,1)', ...
    median(pooledTrialMI,1)', ...
    std(pooledTrialMI,0,1)', ...
    iqr(pooledTrialMI,1)', ...
    min(pooledTrialMI,[],1)', ...
    max(pooledTrialMI,[],1)', ...
    'VariableNames', {'Band','Mean','Median','Std','IQR','Min','Max'});
descFile = fullfile(sessDir, 'descriptives.csv');
    writetable(stats, descFile);

    %% Boxplot
    fig = figure('Visible','off');
    meanMI_perTrial = mean(pooledTrialMI, 2);  % collapse across bands
    boxplot(meanMI_perTrial, tbl.Label, 'Labels',{'Loss','Win'});
    xlabel('Outcome'); ylabel('Mean MI across bands');
    title(sprintf('Win vs Loss Mean MI - %s', lvlName));
    boxPath = fullfile(sessDir,'win_vs_loss_boxplot.png');
    print(fig, boxPath,'-dpng'); close(fig);
    fprintf('Saved boxplot to %s\n', boxPath);

    %% Violin plot
    try
        fig = figure('Visible','off');
        violinplot(pooledTrialMI, tbl.Label);
        set(gca,'XTickLabel',{'Loss','Win'});
        ylabel('MI'); title(sprintf('Win vs Loss MI - %s', lvlName));
        vPath = fullfile(sessDir,'win_vs_loss_violin.png');
        print(fig, vPath,'-dpng'); close(fig);
        fprintf('Saved violin plot to %s\n', vPath);
    catch
        fprintf('Violin plot function not available. Skipping.\n');
    end

    %% Bar chart average MI per band with error bars
    fig = figure('Visible','off');
    bar(1:nBands, mean(pooledTrialMI));
    hold on;
    errorbar(1:nBands, mean(pooledTrialMI), std(pooledTrialMI)/sqrt(size(pooledTrialMI,1)), '.k');
    xlabel('Gamma Band'); ylabel('Mean MI');
    title(sprintf('Average MI per Band - %s', lvlName));
    barPath = fullfile(sessDir, 'avg_mi_bar.png');
    print(fig, barPath,'-dpng'); close(fig);
    fprintf('Saved bar chart to %s\n', barPath);

    %% Summary log
    logPath = fullfile(sessDir,'summary_log.txt');
    fid = fopen(logPath,'w');
    fprintf(fid,'Subject: %s\nSession: %d\nArea: %s\n',subjectID, sessionNum, lvlName);
    fprintf(fid,'Num Trials: %d\n', size(pooledTrialMI,1));
    fprintf(fid,'Win trials: %d\nLoss trials: %d\n', sum(winIdx), sum(lossIdx));
    fprintf(fid,'\nLogistic Regression Coefficients:\n');
    disp(coefTbl, fid);
    fprintf(fid,'\nAUC: %.4f\n', AUC);
    fprintf(fid,'Timestamp: %s\n', datestr(now));
    fclose(fid);
    fprintf('Saved summary log to %s\n', logPath);
end
end
