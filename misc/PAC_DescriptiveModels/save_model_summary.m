function save_model_summary(mdl, resultDir, featureNames)
%SAVE_MODEL_SUMMARY Save logistic regression results and diagnostics.
%   save_model_summary(mdl, resultDir, featureNames) saves coefficient
%   information, diagnostic plots, and model statistics for the provided
%   GeneralizedLinearModel object MDL into RESULTDIR. FEATURENAMES should be
%   a cell array containing the names of the predictor variables in the
%   order used to train MDL.
%
%   The function writes:
%       coefficients.csv  - table of coefficients with SE and p-values
%       coefficients.mat  - MATLAB version of the table
%       glm_model.mat     - the full GeneralizedLinearModel object
%       model_summary.txt - textual output from the model
%       coef_bar.png      - bar chart of coefficients with 95% CI
%       roc_curve.png     - ROC curve for model predictions
%       confusion_matrix.png - Confusion matrix at 0.5 threshold
%       residuals.png     - (optional) residual diagnostic plot
%
%   RESULTDIR must exist prior to calling this function.

if nargin < 3 || isempty(featureNames)
    featureNames = mdl.PredictorNames;
end

% Save coefficient table
coefTbl = mdl.Coefficients;
coefCSV = fullfile(resultDir, 'coefficients.csv');
coefMAT = fullfile(resultDir, 'coefficients.mat');
writecell([coefTbl.Properties.VariableNames; table2cell(coefTbl)], coefCSV);
save(coefMAT, 'coefTbl');

% Save model object
save(fullfile(resultDir, 'glm_model.mat'), 'mdl');

% Create textual summary
summaryStr = evalc('disp(mdl)');
summaryFile = fullfile(resultDir, 'model_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, '%s\n', summaryStr);
% Add goodness-of-fit metrics
pseudoR2 = 1 - mdl.Deviance/mdl.NullDeviance;
fprintf(fid, '\nAIC: %.4f\nDeviance: %.4f\nPseudo-R2: %.4f\n', mdl.ModelCriterion.AIC, mdl.Deviance, pseudoR2);
fclose(fid);

% Plot coefficients with 95% CI
fig = figure('Visible', 'off');
coef = coefTbl.Estimate;
ci = coefTbl.SE * 1.96; % approximate 95% CI
bar(1:numel(coef), coef);
hold on;
errorbar(1:numel(coef), coef, ci, '.k', 'LineWidth', 1);
set(gca, 'XTick', 1:numel(coef), 'XTickLabel', coefTbl.Properties.RowNames, 'XTickLabelRotation', 45);
ylabel('Coefficient Estimate');
title('Logistic Regression Coefficients');
box off;
print(fig, fullfile(resultDir, 'coef_bar.png'), '-dpng');
close(fig);

% ROC curve
fig = figure('Visible', 'off');
[X,Y,~,AUC] = perfcurve(mdl.Variables.(mdl.ResponseName), mdl.Fitted.Probability, 1);
plot(X,Y); grid on;
xlabel('False positive rate'); ylabel('True positive rate');
title(sprintf('ROC Curve (AUC = %.3f)', AUC));
print(fig, fullfile(resultDir, 'roc_curve.png'), '-dpng');
close(fig);

% Confusion matrix at 0.5 threshold
predLabels = mdl.Fitted.Probability >= 0.5;
cm = confusionmat(mdl.Variables.(mdl.ResponseName), predLabels);
fig = figure('Visible', 'off');
confusionchart(cm, {'0','1'});
print(fig, fullfile(resultDir, 'confusion_matrix.png'), '-dpng');
close(fig);

% Residuals diagnostics
fig = figure('Visible', 'off');
plotResiduals(mdl, 'fitted');
print(fig, fullfile(resultDir, 'residuals.png'), '-dpng');
close(fig);
end
