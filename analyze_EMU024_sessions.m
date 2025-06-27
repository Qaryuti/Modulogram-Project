function analyze_EMU024_sessions(baseDir, outDir)
% ANALYZE_EMU024_SESSIONS Analyze modulogram outputs for EMU024 across three sessions.
%   analyze_EMU024_sessions(BASEDIR, OUTDIR) expects BASEDIR to contain folders
%   produced by batch_modulogram.m for EMU024 sessions 1-3 with the win
%   alignment. Each folder must include an AUC_ranking.txt file under
%   EMU024/Session_N and a corresponding EMU024_sessionN_struct.mat file.
%   Results are written to OUTDIR, which defaults to an "analysis" subfolder
%   inside BASEDIR if not specified.
%
%   The script logs progress and errors to analysis_log_sessionN.txt for
%   each session and saves summary figures and correlation tables in the
%   analysis folder.
%
%   Example:
%       analyze_EMU024_sessions('trial-run-3', 'trial-run-3/analysis');

if nargin < 1
    baseDir = 'trial-run-3';
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(baseDir, 'EMU024_analysis');
end

subjectID    = 'EMU024';
sessionNums  = 1:3;
sessionDirs  = arrayfun(@(n) fullfile(baseDir, sprintf('%s_win_mod30-120_session%d_bw5_bins36', subjectID, n)), ...
                        sessionNums, 'UniformOutput', false);
analysisDir  = outDir;
if ~exist(analysisDir, 'dir'); mkdir(analysisDir); end

%% Storage for later comparisons
aucTop  = cell(1, numel(sessionNums));
meanMI  = cell(1, numel(sessionNums));
trialMI = cell(1, numel(sessionNums));
roiAvg  = cell(1, numel(sessionNums));
roiChLabels = cell(1, numel(sessionNums));
freqs   = [];
binCenters = [];
logFIDs = zeros(1, numel(sessionNums));

%% --- Per-session processing ---
for i = 1:numel(sessionNums)
    ses = sessionNums(i);
    sesDir = sessionDirs{i};
    logPath = fullfile(analysisDir, sprintf('analysis_log_session%d.txt', ses));
    logFID = fopen(logPath, 'a');
    logFIDs(i) = logFID;
    logf(logFID, '\n=== Starting analysis for Session %d ===\n', ses);

    structFile = fullfile(sesDir, sprintf('%s_session%d_struct.mat', subjectID, ses));
    if ~exist(structFile, 'file')
        logf(logFID, 'Struct file missing: %s\n', structFile);
        fclose(logFID); logFIDs(i) = -1; continue; end
    try
        data = load(structFile);
    catch ME
        logf(logFID, 'Error loading struct: %s\n', ME.message);
        fclose(logFID); logFIDs(i) = -1; continue;
    end
    try
        sessStruct = data.allData.(subjectID).session(ses).alignment.win.channel;
    catch ME
        logf(logFID, 'Malformed structure: %s\n', ME.message);
        fclose(logFID); logFIDs(i) = -1; continue;
    end
    chNames = fieldnames(sessStruct);
    logf(logFID, 'Found %d channels. Listing regions...\n', numel(chNames));
    anatomy = cell(numel(chNames),1);
    for c = 1:numel(chNames)
        region = sessStruct.(chNames{c}).anatomicalRegion;
        anatomy{c} = region;
        logf(logFID, '  %s - %s\n', chNames{c}, region);
    end

    %% Read AUC ranking
    aucFile = fullfile(sesDir, subjectID, sprintf('Session_%d', ses), 'AUC_ranking.txt');
    if exist(aucFile, 'file')
        try
            T = readtable(aucFile, 'FileType','text', 'ReadVariableNames', false);
            T.Properties.VariableNames = {'Channel','AUC'};
            [~, idx] = sort(T.AUC, 'descend');
            topIdx = idx(1:min(10, height(T)));
            aucTop{i} = T.AUC(topIdx);
            logf(logFID, 'Top 10 channels by AUC:\n');
            for k = 1:numel(topIdx)
                logf(logFID, '  %s = %.4f\n', T.Channel{topIdx(k)}, T.AUC(topIdx(k)));
            end
        catch ME
            logf(logFID, 'Unable to read AUC ranking: %s\n', ME.message);
            aucTop{i} = [];
        end
    else
        logf(logFID, 'AUC ranking not found: %s\n', aucFile);
        aucTop{i} = [];
    end

    %% Extract MI data
    nFreq = numel(sessStruct.(chNames{1}).centerFreqs);
    if isempty(freqs)
        freqs = sessStruct.(chNames{1}).centerFreqs(:);
        binCenters = sessStruct.(chNames{1}).binCenters(:);
    end
    miMat = nan(nFreq, numel(chNames));
    nTrials = size(sessStruct.(chNames{1}).MI_per_trial, 2);
    trialMat = zeros(nFreq, nTrials);
    for c = 1:numel(chNames)
        miMat(:,c) = sessStruct.(chNames{c}).finalAggMI(:);
        trialMat = trialMat + sessStruct.(chNames{c}).MI_per_trial;
    end
    meanMI{i}  = mean(miMat, 2, 'omitnan');
    trialMI{i} = trialMat ./ numel(chNames);

    %% Posterior Cingulate ROI
    roiCh = chNames(strcmp(anatomy, 'Posterior Cingulate'));
    if ~isempty(roiCh)
        nPhase = size(sessStruct.(roiCh{1}).finalAggAmplP, 2);
        roiCube = nan(nFreq, nPhase, numel(roiCh));
        for r = 1:numel(roiCh)
            roiCube(:,:,r) = sessStruct.(roiCh{r}).finalAggAmplP;
        end
        roiAvg{i} = mean(roiCube, 3, 'omitnan');
        roiChLabels{i} = roiCh;
        logf(logFID, 'Posterior Cingulate channels: %s\n', strjoin(roiCh, ', '));
    else
        logf(logFID, 'No Posterior Cingulate channels found.\n');
        roiAvg{i} = [];
        roiChLabels{i} = {};
    end

    % keep logFID open for cross-session summaries
end

%% --- Correlations of top AUCs ---
corFile = fullfile(analysisDir, 'AUC_correlations.txt');
cfid = fopen(corFile, 'w');
aucPairs = {[1 2], [1 3], [2 3]};
for p = 1:numel(aucPairs)
    i = aucPairs{p}(1); j = aucPairs{p}(2);
    if isempty(aucTop{i}) || isempty(aucTop{j})
        rho = NaN;
    else
        rho = corr(aucTop{i}, aucTop{j}, 'Type','Spearman', 'Rows','pairwise');
    end
    line = sprintf('Session %d vs %d Spearman rho: %.4f\n', i, j, rho);
    fprintf(cfid, line);
    if logFIDs(i) > 0
        logf(logFIDs(i), line);
    end
    if logFIDs(j) > 0 && j ~= i
        logf(logFIDs(j), line);
    end
end
fclose(cfid);

%% --- MI curves across sessions ---
fig = figure('Visible','off'); hold on;
cols = lines(numel(meanMI));
anyLines = false;
for i = 1:numel(meanMI)
    if ~isempty(meanMI{i})
        plot(freqs, meanMI{i}, 'Color', cols(i,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('Session %d', sessionNums(i)));
        anyLines = true;
    end
end
xlabel('Gamma Frequency (Hz)');
ylabel('Mean MI');
if anyLines
    legend('Location','best');
end
title('MI Across Sessions');
miCurveFile = fullfile(analysisDir, 'MI_across_sessions.png');
saveas(fig, miCurveFile); close(fig);
for i = 1:numel(logFIDs)
    if logFIDs(i) > 0
        logf(logFIDs(i), 'Saved figure: %s\n', miCurveFile);
    end
end

%% --- Gamma band statistics and boxplot ---
bands = [30 50; 50 70; 70 90; 90 110];
bandMI = zeros(numel(sessionNums), size(bands,1));
for b = 1:size(bands,1)
    idx = freqs >= bands(b,1) & freqs < bands(b,2);
    for i = 1:numel(meanMI)
        if ~isempty(meanMI{i})
            bandMI(i,b) = mean(meanMI{i}(idx), 'omitnan');
        else
            bandMI(i,b) = NaN;
        end
    end
end
bandMeans = mean(bandMI, 1, 'omitnan');
bandStd   = std(bandMI, 0, 1, 'omitnan');

fig = figure('Visible','off');
boxplot(bandMI, 'Labels', {'30-50','50-70','70-90','90-110'});
ylabel('MI');
title('MI per Gamma Band Across Sessions');
boxFile = fullfile(analysisDir, 'MI_band_boxplot.png');
saveas(fig, boxFile); close(fig);

for i = 1:numel(logFIDs)
    if logFIDs(i) > 0
        logf(logFIDs(i), 'Saved figure: %s\n', boxFile);
        logf(logFIDs(i), 'Band statistics:\n');
        for b = 1:size(bands,1)
            logf(logFIDs(i), '  %d-%d Hz mean=%.4f std=%.4f\n', bands(b,1), bands(b,2), bandMeans(b), bandStd(b));
        end
    end
end

%% --- Posterior Cingulate ROI average modulogram ---
validRoi = ~cellfun(@isempty, roiAvg);
if any(validRoi)
    nPhase = size(roiAvg{find(validRoi,1)}, 2);
    roiStack = nan(numel(freqs), nPhase, sum(validRoi));
    idx = 1;
    for i = 1:numel(sessionNums)
        if validRoi(i)
            roiStack(:,:,idx) = roiAvg{i};
            idx = idx + 1;
        end
    end
    roiMean = mean(roiStack, 3, 'omitnan');
    [PH, FR] = meshgrid(binCenters, freqs);
    fig = figure('Visible','off');
    surf(PH, FR, roiMean, 'EdgeColor','none');
    xlabel('\theta phase (rad)'); ylabel('\gamma frequency (Hz)');
    zlabel('Normalized amplitude');
    title('Posterior Cingulate Average Modulogram');
    view([45 25]); colormap turbo; shading interp; colorbar;
    roiFile = fullfile(analysisDir, 'ROI_PosteriorCingulate.png');
    saveas(fig, roiFile); close(fig);
    for i = 1:numel(logFIDs)
        if logFIDs(i) > 0
            logf(logFIDs(i), 'Saved ROI figure: %s\n', roiFile);
            logf(logFIDs(i), 'ROI channels this session: %s\n', strjoin(roiChLabels{i}, ', '));
        end
    end
end

%% --- Trial-by-trial MI drift ---
for i = 1:numel(trialMI)
    if isempty(trialMI{i}) || logFIDs(i) <= 0
        continue; end
    nTrials = size(trialMI{i}, 2);
    meanPerTrial = mean(trialMI{i}, 1, 'omitnan');
    mdl = fitlm((1:nTrials)', meanPerTrial');
    slope = mdl.Coefficients.Estimate(2);
    pval  = mdl.Coefficients.pValue(2);
    note = '';
    if pval < 0.05
        note = ' **SIGNIFICANT**';
    end
    logf(logFIDs(i), 'MI drift slope=%.6f p=%.4f%s\n', slope, pval, note);
end

%% --- k-means clustering on trial-level MI vectors ---
allTrials = [];
for i = 1:numel(trialMI)
    if ~isempty(trialMI{i})
        allTrials = [allTrials, trialMI{i}];
    end
end
allTrials = allTrials'; % trials x freq
if ~isempty(allTrials)
    [idx, C] = kmeans(allTrials, 2, 'Replicates', 5, 'Display', 'off');
    fig = figure('Visible','off');
    plot(freqs, C(1,:), 'LineWidth', 2); hold on;
    plot(freqs, C(2,:), 'LineWidth', 2);
    legend({'Cluster 1','Cluster 2'}); xlabel('Gamma Frequency (Hz)'); ylabel('MI');
    title('MI Cluster Centers');
    clusterFile = fullfile(analysisDir, 'MI_clusters.png');
    saveas(fig, clusterFile); close(fig);
    for i = 1:numel(logFIDs)
        if logFIDs(i) > 0
            logf(logFIDs(i), 'Saved cluster figure: %s\n', clusterFile);
            logf(logFIDs(i), 'Cluster 1 center: %s\n', mat2str(C(1,:),4));
            logf(logFIDs(i), 'Cluster 2 center: %s\n', mat2str(C(2,:),4));
        end
    end
end

%% --- Close logs ---
for i = 1:numel(logFIDs)
    if logFIDs(i) > 0
        logf(logFIDs(i), 'Analysis completed successfully at %s\n', datestr(now));
        fclose(logFIDs(i));
    end
end
end

function logf(fid, fmt, varargin)
    timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
    if fid > 0
        fprintf(fid, ['[' timestamp '] ' fmt], varargin{:});
    end
end
