function outputStruct = create_single_modulogram(config, subjectID, sessionNum, channelLabel, alignment, Fs, filters, trial_times, trial_words, anatomicalRegion)
%CREATE_SINGLE_MODULOGRAM Compute aggregated modulogram for one channel.

outputStruct = [];

% Load channel data
baseDir = fullfile(config.dataDir, '1_formatted');
channelFile = fullfile(baseDir, subjectID, 'separate_channel_files', ...
    sprintf('%s_MAD_SES%d_ch%s.mat', subjectID, sessionNum, channelLabel));
if ~exist(channelFile, 'file')
    fprintf('Missing channel file for %s, Session %d, Channel %s\n', subjectID, sessionNum, channelLabel);
    return;
end
load(channelFile, 'data');
data = (data - mean(data)) / std(data);

% Alignment times
try
    [align_times, ~] = get_align_times(filters, trial_times, trial_words, alignment);
catch ME
    fprintf('Error in get_align_times for %s, Session %d, Channel %s: %s\n', subjectID, sessionNum, channelLabel, ME.message);
    return;
end
align_times = align_times(~isnan(align_times));
if isempty(align_times)
    fprintf('No valid alignment times for %s, Session %d, Channel %s\n', subjectID, sessionNum, channelLabel);
    return;
end
align_times = round(align_times * Fs);
nTrials = numel(align_times);

% Filtering using unified FIR implementation
theta_full = apply_fir_filter(data, Fs, config.modWaveRanges.Theta, config.fir_order);
comb_full  = apply_fir_filter(data, Fs, config.modulatedRange, config.fir_order);

stepSize  = 3;  % Gamma band step size in Hz
bandwidth = config.bandwidth;  % +/- range around each center frequency
centerFreqs   = config.modulatedRange(1):stepSize:config.modulatedRange(2);
numGammaBands = numel(centerFreqs);
subBand_full  = cell(1, numGammaBands);

for iG = 1:numGammaBands
    cfreq = centerFreqs(iG);
    lowBound  = max(cfreq - bandwidth, 1);
    highBound = min(cfreq + bandwidth, Fs/2 - 1);
    subBand_full{iG} = apply_fir_filter(data, Fs, [lowBound, highBound], config.fir_order);
end

% Snip trials
[thetaSnips, ~] = snip_alignment_allTrials(theta_full, align_times, Fs, config.windowParams);
[combSnips, ~] = snip_alignment_allTrials(comb_full, align_times, Fs, config.windowParams);
subBandSnips = cell(numGammaBands,1);
for iG = 1:numGammaBands
    [snips_i, ~] = snip_alignment_allTrials(subBand_full{iG}, align_times, Fs, config.windowParams);
    subBandSnips{iG} = snips_i;
end
if isempty(thetaSnips)
    fprintf('No valid trials after snipping for %s, Session %d, Channel %s\n', subjectID, sessionNum, channelLabel);
    return;
end

% Compute modulogram per trial and average
allAmplP_sub = zeros(numGammaBands, nTrials, config.numBins);
allMI_sub = zeros(numGammaBands, nTrials);
for iG = 1:numGammaBands
    for t = 1:nTrials
        pac_trial.theta_wave = thetaSnips{t};
        pac_trial.constituent_gamma_waves = {subBandSnips{iG}{t}};
        pac_trial.center_freqs = centerFreqs(iG);
        [AmplP_trial, MI_trial, binCenters_trial, freq_array_trial] = compute_modulogram_analysis(pac_trial, Fs, config.numBins);
        allAmplP_sub(iG,t,:) = AmplP_trial;
        allMI_sub(iG,t) = MI_trial;
    end
end
aggAmplP_sub = squeeze(mean(allAmplP_sub,2));
aggMI_sub = mean(allMI_sub,2);
normalized_aggAmplP_sub = zeros(size(aggAmplP_sub));
for iG = 1:size(aggAmplP_sub,1)
    currentRow = aggAmplP_sub(iG,:);
    normalized_aggAmplP_sub(iG,:) = currentRow / mean(currentRow);
end
finalAggAmplP = normalized_aggAmplP_sub;
finalAggMI = aggMI_sub;

resultsDir = fullfile(config.resultsDir, subjectID, sprintf('Session_%d', sessionNum));
if ~exist(resultsDir,'dir'); mkdir(resultsDir); end
outFile = fullfile(resultsDir, sprintf('CH%s_Modulogram.png', channelLabel));
plot_modulogram_figure(finalAggAmplP, finalAggMI, binCenters_trial, centerFreqs, outFile, anatomicalRegion);

outputStruct = struct( ...
    'MI_per_trial', allMI_sub, ...
    'AmplP_per_trial', allAmplP_sub, ...
    'finalAggMI', finalAggMI, ...
    'finalAggAmplP', finalAggAmplP, ...
    'centerFreqs', centerFreqs, ...
    'binCenters', binCenters_trial, ...
    'anatomicalRegion', anatomicalRegion, ...
    'nTrials', nTrials );
end
