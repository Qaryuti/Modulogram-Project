function testt(config)
    % RUN_MODULOGRAM_PIPELINE_WRAPPER
    % This wrapper loops over subjects, sessions, alignments, and channels.
    % It loads the Setup file for each subject/session to obtain the full
    % channel list (elec_ind) and then calls create_single_modulogram for each channel.
    
    % Define subject-specific parameters in a mapping.
    subjectParams = struct(...
        'EMU001', struct('num_sessions', 3, 'Fs', 1000), ...
        'EMU024', struct('num_sessions', 3, 'Fs', 2048), ...
        'EMU025', struct('num_sessions', 2, 'Fs', 2048), ...
        'EMU030', struct('num_sessions', 2, 'Fs', 2048), ...
        'EMU036', struct('num_sessions', 4, 'Fs', 2048), ... % session 2 will be skipped
        'EMU037', struct('num_sessions', 4, 'Fs', 2048), ...
        'EMU038', struct('num_sessions', 1, 'Fs', 2048), ...
        'EMU039', struct('num_sessions', 4, 'Fs', 2048), ...
        'EMU040', struct('num_sessions', 2, 'Fs', 2048), ...
        'EMU041', struct('num_sessions', 9, 'Fs', 2048), ...
        'EMU047', struct('num_sessions', 1, 'Fs', 2048), ...
        'EMU051', struct('num_sessions', 1, 'Fs', 2048) ...
    );
    
    for subjIdx = 1:numel(config.subjects)
        subjectID = config.subjects{subjIdx};
        fprintf('Processing Subject: %s\n', subjectID);
        
        if ~isfield(subjectParams, subjectID)
            fprintf('  No parameters defined for subject %s, skipping.\n', subjectID);
            continue;
        end
        
        num_sessions = subjectParams.(subjectID).num_sessions;
        Fs = subjectParams.(subjectID).Fs;
        
        for sesnum = 1:num_sessions
            fprintf('  Session: %d/%d\n', sesnum, num_sessions);
            
            % Special-case: skip session 2 for EMU036.
            if strcmpi(subjectID, 'EMU036') && sesnum == 2
                fprintf('    Skipping session 2 for %s.\n', subjectID);
                continue;
            end
            
            % Build the Setup file path.
            data_base_dir = fullfile(config.dataDir, '1_formatted');
            if strcmpi(config.reference, 'Ground')
                setupFile = fullfile(data_base_dir, subjectID, sprintf('%s_MAD_SES%d_Setup.mat', subjectID, sesnum));
            else
                setupFile = fullfile(data_base_dir, subjectID, sprintf('%s_MAD_SES%d_Setup_%s.mat', subjectID, sesnum, config.reference));
            end
            
            if ~exist(setupFile, 'file')
                fprintf('    Setup file missing for %s, Session %d\n', subjectID, sesnum);
                continue;
            end
            
            fprintf('    Loading Setup file: %s\n', setupFile);
            load(setupFile, 'filters', 'trial_times', 'trial_words', 'elec_area', 'elec_ind');
            % Now, elec_ind holds the exhaustive list of channels.
            
            % Loop through each alignment in config.
            for alignIdx = 1:numel(config.alignments)
                alignment = config.alignments{alignIdx};
                fprintf('    Alignment: %s\n', alignment);
                
                % Loop through each channel from elec_ind.
                for chIdx = 1:numel(elec_ind)
                    channelLabel = sprintf('%03d', elec_ind(chIdx));
                    
                    % Extract anatomical region for the current channel.
                    if exist('elec_area', 'var') && ~isempty(elec_area)
                        if iscell(elec_area)
                            anatomicalRegion = elec_area{chIdx};
                        else
                            % If elec_area is a matrix (or char array), get the row corresponding to the channel.
                            anatomicalRegion = strtrim(elec_area(chIdx, :));
                        end
                    else
                        anatomicalRegion = 'Unknown';
                    end

                    % Call create_single_modulogram with the additional anatomicalRegion parameter.
                    create_single_modulogram(config, subjectID, sesnum, channelLabel, alignment, Fs, filters, trial_times, trial_words, anatomicalRegion);
                end
            end
        end
    end
end

function create_single_modulogram(config, subjectID, sessionNum, channelLabel, alignment, Fs, filters, trial_times, trial_words, anatomicalRegion)
    % CREATE_SINGLE_MODULOGRAM computes an aggregate modulogram for a given channel
    % using sub-band PAC analysis.
    %
    % Inputs:
    %   config       - Configuration struct with fields:
    %                  .dataDir        : Base directory for data.
    %                  .resultsDir     : Directory where results (plots) are saved.
    %                  .modulatedRange : Frequency range for gamma (e.g., [30 130]).
    %                  .modWaveRanges  : Structure with .Theta (e.g., [4 8]).
    %                  .bandwidth      : Bandwidth for subdividing gamma (e.g., 5).
    %                  .numBins        : Number of phase bins for modulogram (e.g., 64).
    %                  .windowParams   : Struct with window parameters (e.g., window: [0.5, 0.5]).
    %   subjectID    - e.g., 'EMU024'
    %   sessionNum   - Session number (integer)
    %   channelLabel - A string label for the channel (e.g., '003' for Ground)
    %   alignment    - Alignment condition (e.g., 'loss')
    %   Fs           - Sampling rate for this subject/session.
    %   filters, trial_times, trial_words - Variables loaded from the Setup file.
    %
    % This function:
    %   1. Loads and normalizes the channel data.
    %   2. Extracts alignment times.
    %   3. Filters the entire time series for theta (phase), combined gamma,
    %      and gamma sub-bands.
    %   4. Snips the filtered signals around each alignment event.
    %   5. Computes the modulogram analysis for each gamma sub-band trial-by-trial,
    %      averages across trials for each sub-band, and then averages across sub-bands.
    %   6. Plots the aggregated modulogram (heatmap and MI line plot) and saves the plot.
    
    %% 1. Load and Normalize Data
    data_base_dir = fullfile(config.dataDir, '1_formatted');
    channelFile = fullfile(data_base_dir, subjectID, 'separate_channel_files', ...
        sprintf('%s_MAD_SES%d_ch%s.mat', subjectID, sessionNum, channelLabel));
    if ~exist(channelFile, 'file')
        fprintf('Missing channel file for %s, Session %d, Channel %s\n', subjectID, sessionNum, channelLabel);
        return;
    end
    load(channelFile, 'data');  % Assumes variable 'data' exists
    data = (data - mean(data)) / std(data);  % Normalize (z-score)
    
    %% 2. Extract Alignment Times
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
    fprintf('Aggregating %d trials for %s, Session %d, Channel %s\n', nTrials, subjectID, sessionNum, channelLabel);
    
    %% 3. Filter Entire Time Series
    % (a) Theta wave (for phase extraction).
    theta_full = butterworth_filter(data, Fs, config.modWaveRanges.Theta);
    % (b) Combined gamma wave.
    comb_full = butterworth_filter(data, Fs, config.modulatedRange);
    
    % (c) Gamma sub-band signals.
    stepSize = 2 * config.bandwidth;
    centerFreqs = config.modulatedRange(1):stepSize:config.modulatedRange(2);
    numGammaBands = numel(centerFreqs);
    subBand_full = cell(1, numGammaBands);
    for iG = 1:numGammaBands
        cfreq = centerFreqs(iG);
        lowBound = cfreq - config.bandwidth;
        highBound = cfreq + config.bandwidth;
        subBand_full{iG} = butterworth_filter(data, Fs, [lowBound, highBound]);
    end

    %% Debug Print Statements After Filtering

    fprintf('--- Filtering Summary ---\n');
    fprintf('Theta wave: Length = %d samples, first 5 values: %s\n', length(theta_full), mat2str(theta_full(1:min(5, length(theta_full))), 4));
    fprintf('Combined gamma wave: Length = %d samples, first 5 values: %s\n', length(comb_full), mat2str(comb_full(1:min(5, length(comb_full))), 4));

    fprintf('Gamma sub-band signals: %d sub-bands\n', numGammaBands);
    fprintf('Center frequencies: %s\n', mat2str(centerFreqs));

    for iG = 1:numGammaBands
        currentBand = subBand_full{iG};
        fprintf('Sub-band %d (Center frequency = %d Hz): Length = %d samples, first 5 values: %s\n', ...
            iG, centerFreqs(iG), length(currentBand), mat2str(currentBand(1:min(5, length(currentBand))), 4));
    end
    fprintf('--- End of Filtering Summary ---\n');

    
    %% 4. Snip Filtered Signals
    % Here we assume that snip_alignment_allTrials returns a cell array of trials.
    [thetaSnips, ~] = snip_alignment_allTrials(theta_full, align_times, Fs, config.windowParams);
    [combSnips, ~]  = snip_alignment_allTrials(comb_full, align_times, Fs, config.windowParams);
    subBandSnips = cell(numGammaBands, 1);
    for iG = 1:numGammaBands
        [snips_i, ~] = snip_alignment_allTrials(subBand_full{iG}, align_times, Fs, config.windowParams);
        subBandSnips{iG} = snips_i;
    end
    
    if isempty(thetaSnips)
        fprintf('No valid trials after snipping for %s, Session %d, Channel %s\n', subjectID, sessionNum, channelLabel);
        return;
    end

    %% Debug: Print Snipping Summary
    fprintf('--- Snipping Summary ---\n');
    if ~isempty(thetaSnips)
        fprintf('Theta Snippets: %d valid trials extracted. Each snippet length: %d samples.\n', ...
            numel(thetaSnips), length(thetaSnips{1}));
    else
        fprintf('No valid theta snippets extracted.\n');
    end

    if ~isempty(combSnips)
        fprintf('Combined Gamma Snippets: %d valid trials extracted. Each snippet length: %d samples.\n', ...
            numel(combSnips), length(combSnips{1}));
    else
        fprintf('No valid combined gamma snippets extracted.\n');
    end

    for iG = 1:numGammaBands
        if ~isempty(subBandSnips{iG})
            fprintf('Sub-band %d (Center = %d Hz): %d valid trials extracted, snippet length: %d samples.\n', ...
                iG, centerFreqs(iG), numel(subBandSnips{iG}), length(subBandSnips{iG}{1}));
        else
            fprintf('Sub-band %d (Center = %d Hz): No valid trials extracted.\n', iG, centerFreqs(iG));
        end
    end
    fprintf('--- End of Snipping Summary ---\n');

    
   %% 5. Compute Modulation Analysis per Sub-Band and Average Across Trials
% We process each gamma sub-band separately.
allAmplP_sub = zeros(numGammaBands, nTrials, config.numBins); % [subBand x trial x numBins]
allMI_sub = zeros(numGammaBands, nTrials);                     % [subBand x trial]

for iG = 1:numGammaBands
    for t = 1:nTrials
        % Build a temporary PAC signals struct for trial t and sub-band iG.
        pac_trial.theta_wave = thetaSnips{t};
        pac_trial.constituent_gamma_waves = { subBandSnips{iG}{t} };
        pac_trial.center_freqs = centerFreqs(iG);  % scalar for this sub-band
        % Compute modulogram analysis for this trial.
        [AmplP_trial, MI_trial, binCenters_trial, freq_array_trial] = compute_modulogram_analysis(pac_trial, Fs, config.numBins);
        allAmplP_sub(iG, t, :) = AmplP_trial; % AmplP_trial: [1 x numBins]
        allMI_sub(iG, t) = MI_trial;            % MI_trial: scalar
    end
end

% Average the modulogram outputs across trials for each sub-band.
aggAmplP_sub = squeeze(mean(allAmplP_sub, 2));  % [numGammaBands x numBins]
aggMI_sub = mean(allMI_sub, 2);                   % [numGammaBands x 1]

% --- Step 1 (revised): normalize each sub-band so its mean is 1 ---
normAmp       = zeros(size(aggAmplP_sub));
for iG = 1:size(aggAmplP_sub,1)
  normAmp(iG,:) = aggAmplP_sub(iG,:) / mean(aggAmplP_sub(iG,:));
end
finalAggAmplP = normAmp;    % heatmap uses mean‑normalized amplitudes



% ——————————————————————————————————————
% Subritzky‑Katz et al. first-order sinusoid fit per band :contentReference[oaicite:2]{index=2}&#8203;:contentReference[oaicite:3]{index=3}
nbins   = numel(binCenters_trial);             % # of phase bins
nBands  = size(P_sub,1);                       % # of gamma sub-bands

% Pre-build design matrix [a0 + a1*cos(ϕ) + b1*sin(ϕ)]
X = [ ones(nbins,1), cos(binCenters_trial(:)), sin(binCenters_trial(:)) ];  % [nbins×3]

rho_sin = zeros(nBands,1);
depth   = zeros(nBands,1);

for b = 1:nBands
  A      = P_sub(b,:)';                        % [nbins×1] distribution for band b
  coeffs = X \ A;                              % solve [a0;a1;b1]
  fit   = X * coeffs;                          % fitted sinusoid at each ϕ

  % eq (2): fit error
  rho_sin(b) = sqrt( mean((A - fit).^2) );

  % modulation depth (peak‑to‑trough)
  depth(b)   = max(fit) - min(fit);
end

% Use "depth" (or 1 - rho_sin, or rho_sin itself) as your per‑band MI:
MI_perBand = depth;  
% ——————————————————————————————————————



% --- Step 2: fit a first‑order sinusoid to each row and pull out depth & error ---
nBands    = size(P_sub,1);
binCenters = binCenters_trial;            % 1×numBins vector from compute_modulogram_analysis
nbins     = numel(binCenters);

depth     = zeros(nBands,1);
rho_sin   = zeros(nBands,1);

% Build once: design matrix for [a0 + a1 cosϕ + b1 sinϕ]
X = [ ones(nbins,1), cos(binCenters(:)), sin(binCenters(:)) ];  % nbins×3

for iG = 1:nBands
  Amp      = P_sub(iG,:)';               % nbins×1
  coeffs   = X \ Amp;                    % [a0; a1; b1]
  fitVals  = X * coeffs;                 % nbins×1 fitted sinusoid
  rho_sin(iG) = sqrt( mean((Amp - fitVals).^2) );
  depth(iG)   = max(fitVals) - min(fitVals);
end


    
    %% 6. Build Final PAC Signals Structure (Optional)
    % You might store the aggregated values along with centerFreqs.
    pac_signals_agg = struct();
    pac_signals_agg.theta_wave = mean(cell2mat(thetaSnips), 2);  % aggregate of raw theta (if desired)
    pac_signals_agg.combined_gamma_wave = mean(cell2mat(combSnips), 2); % not used further here
    pac_signals_agg.constituent_gamma_waves = aggAmplP_sub; % aggregated per sub-band (optional)
    pac_signals_agg.center_freqs = centerFreqs;
    
%% 7. Plot the Aggregate Modulogram
resultsDir = fullfile(config.resultsDir, subjectID, sprintf('Session_%d', sessionNum));
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
outFile = fullfile(resultsDir, sprintf('CH%s_Modulogram.png', channelLabel));

% Expansive print statements to verify variables.
fprintf('--------------------\n');
fprintf('Aggregated Results for %s, Session %d, Channel %s:\n', subjectID, sessionNum, channelLabel);
fprintf('Number of valid trials: %d\n', nTrials);
fprintf('Number of gamma sub-bands: %d\n', numGammaBands);
fprintf('Center Frequencies: %s\n', mat2str(centerFreqs));
fprintf('Size of final aggregated AmplP (should be [1 x numBins]): %s\n', mat2str(size(finalAggAmplP)));
fprintf('Sinusoidal‑fit depth (per sub-band): %s\n', mat2str(MI_perBand,4));
fprintf('Bin Centers (phase bins) size: %s\n', mat2str(size(binCenters_trial)));
fprintf('Bin Centers: %s\n', mat2str(binCenters_trial, 4));
fprintf('Output file will be saved as: %s\n', outFile);
fprintf('--------------------\n');


fig = figure('Visible','off','Position',[100,100,800,600]);
subplot(2,1,1);
  imagesc(binCenters_trial, centerFreqs, finalAggAmplP);
  axis xy; colorbar;
  xlabel('Theta Phase (rad)');
  ylabel('Gamma Freq (Hz)');
  title(sprintf('Modulogram Heatmap: %s', anatomicalRegion));
subplot(2,1,2);
plot(centerFreqs, MI_perBand, '-o','LineWidth',2);
ylabel('Modulation depth');

  ylim([0 1]);
  xlabel('Gamma Freq (Hz)');
  ylabel('Modulation Depth');
  title('Sinusoidal PAC Depth per Sub-band');
  grid on;
print(fig,'-dpng',outFile);
close(fig);
end

function trialCount = load_data_and_get_trials(config, subjectID, sessionNum, channelLabel, alignment, Fs, filters, trial_times, trial_words)
    % LOAD_DATA_AND_GET_TRIALS loads channel data and uses get_align_times
    % to determine the number of valid trials.
    %
    % It constructs the channel file path based on subject, session, channelLabel,
    % and the reference scheme defined in config, loads the channel data, and then
    % calls get_align_times. It returns the number of valid trial events.
    
    data_base_dir = fullfile(config.dataDir, '1_formatted');
    
    if strcmpi(config.reference, 'Ground')
        channelFile = fullfile(data_base_dir, subjectID, 'separate_channel_files', ...
            sprintf('%s_MAD_SES%d_ch%s.mat', subjectID, sessionNum, channelLabel));
    else
        channelFile = fullfile(data_base_dir, subjectID, 'separate_channel_files', subjectID, ...
            sprintf('%s_MAD_SES%d_ch%s.mat', subjectID, sessionNum, channelLabel));
    end
    
    if ~exist(channelFile, 'file')
        trialCount = 0;
        return;
    end
    
    load(channelFile, 'data');  % Assumes the variable 'data' exists in the file.
    
    try
        [align_times, ~] = get_align_times(filters, trial_times, trial_words, alignment);
    catch ME
        trialCount = 0;
        return;
    end
    
    align_times = align_times(~isnan(align_times));
    align_times = round(align_times * Fs);
    
    trialCount = numel(align_times);
end

function [AmplP_matrix, MI_vector, binCenters, freq_array] = compute_modulogram_analysis(pac_signals, Fs, num_bins)
    % COMPUTE_MODULOGRAM_ANALYSIS
    %
    % Takes the pac_signals struct (from Section 7) and computes:
    %   - Amplitude-phase distribution (AmplP_matrix) across frequencies & phase bins
    %   - A single MI value (MI_vector) per frequency
    %   - Phase bin centers (binCenters)
    %   - The array of gamma frequencies (freq_array)
    %
    % Inputs:
    %   pac_signals.theta_wave            - the low-frequency wave for phase
    %   pac_signals.constituent_gamma_waves - cell array of gamma waves
    %   pac_signals.center_freqs          - array of gamma center frequencies
    %   Fs                                - sampling rate (could also be in pac_signals)
    %   num_bins                          - # of phase bins (e.g. 16)
    %
    % Outputs:
    %   AmplP_matrix - [nFreq x num_bins], amplitude distribution for each freq & bin
    %   MI_vector    - [nFreq x 1], a single modulation index per freq
    %   binCenters   - [1 x num_bins], center of phase bins
    %   freq_array   - the center frequencies from pac_signals

        % Extract fields
        theta_wave   = pac_signals.theta_wave;
        gamma_waves  = pac_signals.constituent_gamma_waves;  % cell array
        freq_array   = pac_signals.center_freqs;

        % 1) Compute theta phase
        theta_phase = angle(hilbert(theta_wave));

        % 2) Create phase bins
        binEdges   = linspace(-pi, pi, num_bins+1);
        binCenters = 0.5*(binEdges(1:end-1) + binEdges(2:end));

        nFreq = numel(gamma_waves);
        AmplP_matrix = zeros(nFreq, num_bins);
        MI_vector    = zeros(nFreq, 1);

        % 3) For each gamma wave, compute amplitude-phase distribution + MI
        for iF = 1:nFreq
            this_gamma = gamma_waves{iF};
            gamma_ampl = abs(hilbert(this_gamma));
        
            % Bin amplitude by theta phase
            amplP = zeros(1, num_bins);
            for bIdx = 1:num_bins
                mask = (theta_phase >= binEdges(bIdx)) & (theta_phase < binEdges(bIdx+1));
                amplP(bIdx) = mean(gamma_ampl(mask));
            end
            AmplP_matrix(iF,:) = amplP;
        
            % Normalize the amplitude distribution to sum to 1 (probability distribution)
            if sum(amplP) > 0
                p = amplP / sum(amplP);
            else
                p = amplP;  % Handle the case where the sum is zero to avoid division by zero
            end
        
            % Compute KL divergence (raw MI)
            uniformDist = ones(1, num_bins) / num_bins;
            MI_raw = sum(p .* log((p + eps) ./ (uniformDist + eps)));
        
            % Optionally, normalize MI so that it ranges between 0 and 1.
            MI_vector(iF) = MI_raw / log(num_bins);
        end        
end

function filtered_data = butterworth_filter(rawData, Fs, freqRange)
    % BUTTERWORTH_FILTER
    % Applies a bandpass Butterworth filter to the input data.
    %
    % INPUTS:
    %   - rawData:     [samples x 1] The time-series data to be filtered.
    %   - Fs:          Sampling frequency (Hz).
    %   - freqRange:   [1 x 2] Frequency range for the bandpass filter ([lowFreq, highFreq]).
    %
    % OUTPUT:
    %   - filtered_data: The bandpass-filtered time-series data.
    
        %% Validate Inputs
        if length(freqRange) ~= 2
            error('freqRange must be a 2-element vector specifying [lowFreq, highFreq].');
        end
        
        lowFreq = freqRange(1);
        highFreq = freqRange(2);
    
        if lowFreq <= 0 || highFreq >= Fs / 2
            error('Frequency range must be within (0, Nyquist Frequency).');
        end
    
        %% Design Butterworth Bandpass Filter
        % Normalize frequency range to Nyquist frequency (Fs/2)
        nyquist = Fs / 2;
        normalized_range = [lowFreq, highFreq] / nyquist;
    
        % Design a 2nd-order Butterworth filter
        filter_order = 2;
        [b, a] = butter(filter_order, normalized_range, 'bandpass');
    
        %% Apply the Filter
        % Use filtfilt to apply the filter forward and backward, avoiding phase distortion
        filtered_data = filtfilt(b, a, rawData);
    
end

function [snips, valid_idxs] = snip_alignment_allTrials(data, align_times, Fs, someParams)
    % SNIP_ALIGNMENT_ALLTRIALS wraps snip_alignment and returns a cell array of trials.
    %
    % It calls snip_alignment to extract segments (as a matrix) and then converts each
    % valid column into a separate cell.
    
    % Call your existing snip_alignment function.
    snipped_matrix = snip_alignment(data, align_times, Fs, someParams);
    
    % Determine valid trial columns (columns not all NaN)
    valid_idxs = find(~all(isnan(snipped_matrix), 1));
    
    % Convert valid columns to a cell array.
    snips = cell(1, numel(valid_idxs));
    for i = 1:numel(valid_idxs)
        snips{i} = snipped_matrix(:, valid_idxs(i));
    end
end

function snipped_data = snip_alignment(data, align_times, Fs, someParams)
    % SNIP_ALIGNMENT Extracts data segments based on alignment times and a specified window.
    %
    % INPUTS:
    %   - data: [samples x 1] Time-series data for a single channel.
    %   - align_times: [1 x N] Vector of alignment times in samples.
    %   - Fs: Sampling frequency (Hz).
    %   - someParams: Struct containing the following fields:
    %       * window: [1 x 2] Pre- and post-alignment window in seconds ([pre_time, post_time]).
    %       * excludeNaNTrials (optional): Logical flag to exclude trials with out-of-bounds indices (default: false).
    %
    % OUTPUT:
    %   - snipped_data: [samples x N] Matrix of extracted segments, where N is the number of align_times.
    
        %% Default Parameters
        if ~isfield(someParams, 'window')
            error('The field "window" is required in someParams.');
        end
    
        if ~isfield(someParams, 'excludeNaNTrials')
            someParams.excludeNaNTrials = true; % Default: include NaN trials
        end
    
        %% Compute Window in Samples
        pre_samples  = round(someParams.window(1) * Fs); % Samples before alignment
        post_samples = round(someParams.window(2) * Fs); % Samples after alignment
        segment_length = pre_samples + post_samples + 1; % Total length of each segment
    
        %% Initialize Output Matrix
        num_alignments = length(align_times);
        snipped_data = NaN(segment_length, num_alignments);
    
        %% Extract Segments
        for i = 1:num_alignments
            % Define the segment range around the current alignment
            start_idx = align_times(i) - pre_samples;
            end_idx = align_times(i) + post_samples;
    
            % Check if the segment is within bounds
            if start_idx < 1 || end_idx > length(data)
                if someParams.excludeNaNTrials
                    continue; % Skip out-of-bounds trials
                end
            else
                % Extract the segment
                snipped_data(:, i) = data(start_idx:end_idx);
            end
        end
    
        %% Remove NaN Trials if Requested
        if someParams.excludeNaNTrials
            snipped_data(:, all(isnan(snipped_data), 1)) = []; % Remove columns with all NaN
        end
end 
    
function plot_modulogram_figure(AmplP_matrix, MI_vector, binCenters, freq_array, outFile, anatomicalRegion)
    % Create a new figure with two subplots.
    fig = figure('Visible','off','Position',[100,100,800,600]);

    %anatomicalRegion = "hello";
    
    % Subplot 1: Heatmap of amplitude-phase distribution.
    subplot(2,1,1);
    imagesc('XData', binCenters, 'YData', freq_array, 'CData', AmplP_matrix);
    axis xy;  % Ensure lower frequencies are at the bottom.
    colorbar;
    xlabel('Theta Phase (radians)');
    ylabel('Gamma Center Frequency (Hz)');
    title(sprintf('Modulogram\nElectrode Region: %s', anatomicalRegion));
    
    % Subplot 2: Line plot for the modulation index.
    subplot(2,1,2);
    plot(freq_array, MI_vector, '-o', 'LineWidth', 2);
    xlabel('Gamma Center Frequency (Hz)');
    ylabel('Modulation Index (MI)');
    title('Modulation Index across Gamma Sub-bands');
    grid on;
    
    % Save the figure to a file.
    print(fig, '-dpng', outFile);
    close(fig);
end


