% RUN_WINLOSS_LOGISTIC_WRAPPER
% Execute logistic analysis for EMU024 win/loss data.
%
% This script loops over sessions, merges win and loss modulograms
% for each session, and then runs analyze_winloss_logistic_pipeline.
%
% Requires:
%   - analyze_winloss_logistic_pipeline.m
%   - convert_modulogram_merge.m
% in the same folder.

baseDir = fullfile('/home/alabwam1/Desktop/modulogram_V7/trial-run-3');
resultsDir = fullfile('/home/alabwam1/Desktop/modulogram_V7/analysis/results');

% make Data subfolder if missing
analysisDataDir = fullfile('analysis','Data');
if ~exist(analysisDataDir,'dir')
    mkdir(analysisDataDir);
end

% run for each session
for sesNum = 1:3
    winFolder = sprintf('EMU024_win_mod30-120_session%d_bw5_bins36', sesNum);
    lossFolder = sprintf('EMU024_loss_mod30-120_session%d_bw5_bins36', sesNum);
    
    winMatFile = fullfile(baseDir, winFolder, sprintf('EMU024_session%d_struct.mat', sesNum));
    lossMatFile = fullfile(baseDir, lossFolder, sprintf('EMU024_session%d_struct.mat', sesNum));
    
    fprintf('\n=== Merging session %d ===\n', sesNum);
    sessStructFile = convert_modulogram_merge(winMatFile, lossMatFile);
    
    % run logistic
    analyze_winloss_logistic_pipeline(sessStructFile, resultsDir);
end
