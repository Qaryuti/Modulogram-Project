function logreg_with_phase_bin(modulogramStruct, outcomeLabels)
%% LOGREG_WITH_PHASE_BIN Logistic regression using MI and peak phase bins.
%  LOGREG_WITH_PHASE_BIN(modStruct, labels) extends the basic model by
%  including the peak phase bin index of the amplitude distribution for each
%  frequency band in addition to the modulation index (MI). This allows
%  assessing whether the preferred phase provides extra predictive power.
%
%  Results are written to a timestamped folder under
%  PAC_DescriptiveModels/Results/logreg_phase_<timestamp>/.
%
%  The folder contains coefficient tables, diagnostic plots, ROC curve and
%  confusion matrix, plus a grouped bar plot comparing MI and phase-bin
%  coefficient magnitudes.
%
%  Example:
%      logreg_with_phase_bin(modStruct, labels)
%
%  See also SAVE_MODEL_SUMMARY.

% Validate inputs
MI = modulogramStruct.MI_per_trial;
Amp = modulogramStruct.AmplP_per_trial; % [nBands x nTrials x nPhaseBins]
[nBands, nTrials, nPhaseBins] = size(Amp);
assert(isequal(size(MI,1), nBands) && size(MI,2) == nTrials, 'MI_per_trial size mismatch');
assert(numel(outcomeLabels) == nTrials, 'Outcome labels must match trial count');

% Determine peak phase bin per trial
[~, peakBin] = max(Amp, [], 3);  % [nBands x nTrials]

% Build predictor table
featureNames = [strcat('MI_B', string(1:nBands)), strcat('PhaseBin_B', string(1:nBands))];
X = array2table([MI'; peakBin']', 'VariableNames', featureNames);
T = [X, table(outcomeLabels(:), 'VariableNames', {'Outcome'})];

% Construct formula including all predictors linearly
formula = 'Outcome ~ 1';
for i = 1:nBands
    formula = [formula, ' + MI_B', num2str(i), ' + PhaseBin_B', num2str(i)]; %#ok<AGROW>
end

mdl = fitglm(T, formula, 'Distribution', 'binomial');

% Prepare result directory
timestamp = datestr(now,'yyyy-mm-dd_HH-MM');
resultDir = fullfile('PAC_DescriptiveModels','Results', ['logreg_phase_' timestamp]);
if ~exist(resultDir,'dir'); mkdir(resultDir); end

% Save standard outputs
save_model_summary(mdl, resultDir, featureNames);

% Additional bar plot comparing MI vs Phase coefficients
coefTbl = mdl.Coefficients;
miIdx = contains(coefTbl.Properties.RowNames, 'MI_B');
phaseIdx = contains(coefTbl.Properties.RowNames, 'PhaseBin_B');
fig = figure('Visible', 'off');
bar(categorical({'MI','Phase'}), [mean(abs(coefTbl.Estimate(miIdx))) mean(abs(coefTbl.Estimate(phaseIdx)))]);
ylabel('Mean |Coefficient|');
title('Contribution of MI vs Phase Bin');
print(fig, fullfile(resultDir,'MI_vs_Phase_contribution.png'), '-dpng');
close(fig);
end
