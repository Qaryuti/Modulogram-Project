function logreg_with_interactions(modulogramStruct, outcomeLabels)
%% LOGREG_WITH_INTERACTIONS Logistic regression with MI interactions.
%  LOGREG_WITH_INTERACTIONS(modStruct, labels) fits a logistic regression
%  using modulation index (MI) from each band and also includes pairwise
%  interaction terms between bands (i.e., the product of two band MIs).
%
%  The intent is to evaluate whether co-modulation of frequency bands better
%  predicts the outcome than individual bands alone.
%
%  Results are stored in a timestamped folder under
%  PAC_DescriptiveModels/Results/logreg_interactions_<timestamp>/.
%
%  The folder includes coefficient tables, diagnostic plots, ROC curve,
%  confusion matrix, and the fitted GLM object.
%
%  Example:
%      logreg_with_interactions(modStruct, labels)
%
%  See also SAVE_MODEL_SUMMARY.

% Validate inputs
MI = modulogramStruct.MI_per_trial;
[nBands, nTrials] = size(MI);
assert(numel(outcomeLabels) == nTrials, 'Outcome labels must match trial count');

% Base predictors
featureNames = strcat('MI_B', string(1:nBands));
X = MI'; % nTrials x nBands

% Generate interaction terms
interNames = {};
interData = [];
for i = 1:nBands-1
    for j = i+1:nBands
        interNames{end+1} = sprintf('MI_B%d_x_B%d', i, j); %#ok<AGROW>
        interData = [interData, (MI(i,:).*MI(j,:))']; %#ok<AGROW>
    end
end

Xtbl = array2table([X interData], 'VariableNames', [featureNames, interNames]);
T = [Xtbl, table(outcomeLabels(:), 'VariableNames', {'Outcome'})];

% Construct formula
formula = 'Outcome ~ 1';
for i = 1:numel(featureNames)
    formula = [formula, ' + ', featureNames{i}]; %#ok<AGROW>
end
for i = 1:numel(interNames)
    formula = [formula, ' + ', interNames{i}]; %#ok<AGROW>
end

mdl = fitglm(T, formula, 'Distribution', 'binomial');

% Result directory
timestamp = datestr(now,'yyyy-mm-dd_HH-MM');
resultDir = fullfile('PAC_DescriptiveModels','Results', ['logreg_interactions_' timestamp]);
if ~exist(resultDir,'dir'); mkdir(resultDir); end

% Save outputs
save_model_summary(mdl, resultDir, [featureNames, interNames]);
end
