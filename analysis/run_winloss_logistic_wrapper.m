%RUN_WINLOSS_LOGISTIC_WRAPPER Execute logistic analysis for EMU024 win/loss data.
%
%   This script iterates over the modulogram result folders for
%   EMU024_win_mod30-120_bw5_bins36 and EMU024_loss_mod30-120_bw5_bins36
%   (sessions 1-3). It calls analyze_winloss_logistic_pipeline for each
%   session structure and writes results under analysis/results/.
%
%   Example:
%       run_winloss_logistic_wrapper
%
%   Requires analyze_winloss_logistic_pipeline.m in the same folder.

baseDir = fullfile('trial-run-3');
resultsDir = fullfile('analysis','results');

folders = {
    'EMU024_win_mod30-120_session1_bw5_bins36',
    'EMU024_win_mod30-120_session2_bw5_bins36',
    'EMU024_win_mod30-120_session3_bw5_bins36',
    'EMU024_loss_mod30-120_session1_bw5_bins36',
    'EMU024_loss_mod30-120_session2_bw5_bins36',
    'EMU024_loss_mod30-120_session3_bw5_bins36'
    };

for i = 1:numel(folders)
    tokens = regexp(folders{i}, 'session(\d+)', 'tokens', 'once');
    sessNum = tokens{1};
    matFile = fullfile(baseDir, folders{i}, sprintf('EMU024_session%s_struct.mat', sessNum));
    fprintf('Processing %s\n', matFile);
    analyze_winloss_logistic_pipeline(matFile, resultsDir);
end
