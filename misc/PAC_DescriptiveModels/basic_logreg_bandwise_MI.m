function basic_logreg_bandwise_MI(modulogramStruct, outcomeLabels)
%% BASIC_LOGREG_BANDWISE_MI  Logistic regression on band-wise MI.
%  BASIC_LOGREG_BANDWISE_MI(modStruct, labels) fits a logistic regression
%  predicting LABELS (1=win,0=loss) from the modulation index (MI) in
%  MODSTRUCT.MI_per_trial. Each frequency band contributes one predictor.
%
%  The script validates input dimensions, fits the model, and stores results
%  in a timestamped subdirectory of PAC_DescriptiveModels/Results.
%  Outputs include coefficient tables, diagnostic plots (ROC curve,
%  confusion matrix, residuals), and the fitted GLM object.
%
%  Result files are suitable for quick descriptive analysis of PAC effects.
%
%  Example:
%      basic_logreg_bandwise_MI(allData.EMU024.session(1).alignment.loss.channel.CH057, labels)
%
%  See also SAVE_MODEL_SUMMARY.

% Validate inputs
MI = modulogramStruct.MI_per_trial;
[nBands, nTrials] = size(MI);
assert(numel(outcomeLabels) == nTrials, 'Outcome labels must match trial count');

% Build predictor table
featureNames = strcat('MI_B', string(1:nBands));
X = array2table(MI', 'VariableNames', featureNames);
T = [X, table(outcomeLabels(:), 'VariableNames', {'Outcome'})];

% Fit logistic regression
% Expand formula to include all bands
formula = ['Outcome ~ 1'];
for i = 1:nBands
    formula = [formula, ' + MI_B', num2str(i)]; %#ok<AGROW>
end
mdl = fitglm(T, formula, 'Distribution','binomial');

% Prepare result directory
timestamp = datestr(now,'yyyy-mm-dd_HH-MM');
resultDir = fullfile('PAC_DescriptiveModels','Results', ['basic_MI_' timestamp]);
if ~exist(resultDir,'dir'); mkdir(resultDir); end

% Save outputs
save_model_summary(mdl, resultDir, featureNames);
end
