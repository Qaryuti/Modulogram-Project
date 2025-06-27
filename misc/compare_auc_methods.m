function compare_auc_methods(pathA, pathB)
    % pathA: path to AUC_ranking.txt for method A (e.g., Butterworth)
    % pathB: path to AUC_ranking.txt for method B (e.g., FIR overlapping)

    T1 = readtable(pathA, 'Delimiter', '\t', 'ReadVariableNames', false);
    T2 = readtable(pathB, 'Delimiter', '\t', 'ReadVariableNames', false);

    % Standardize variable names
    T1.Properties.VariableNames = {'Channel', 'AUC_A'};
    T2.Properties.VariableNames = {'Channel', 'AUC_B'};

    % Merge on Channel
    M = innerjoin(T1, T2, 'Keys', 'Channel');

    % Compute Delta and Rank
    M.Delta = M.AUC_B - M.AUC_A;
    M.Rank_A = tiedrank(-M.AUC_A);
    M.Rank_B = tiedrank(-M.AUC_B);

    % Spearman correlation
    rho = corr(M.Rank_A, M.Rank_B, 'type', 'Spearman');
    fprintf('Spearman rank correlation: %.3f\n', rho);
    fprintf('Mean ΔAUC: %.4f\n', mean(abs(M.Delta)));

    % Scatter plot of AUCs
    figure;
    subplot(1,2,1);
    scatter(M.AUC_A, M.AUC_B, 'filled');
    hold on; plot([0 1], [0 1], 'r--');
    xlabel('AUC Method A (Butterworth)');
    ylabel('AUC Method B (FIR)');
    title(sprintf('AUC Comparison (\\rho = %.2f)', rho));
    grid on;

    % Histogram of Delta AUC
    subplot(1,2,2);
    histogram(M.Delta, 15);
    xlabel('ΔAUC (FIR - Butter)');
    ylabel('Count');
    title('Difference in AUC per Channel');
    xline(0, 'r--');

    % Optional: save table
    writetable(M, 'AUC_Comparison_Table.csv');

    print(gcf, 'AUC_Comparison_Figure', '-dpng', '-r300');
end
