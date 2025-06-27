function run_pac_stats(config)
    % Load data
    lossData = load(config.lossPath);
    winData  = load(config.winPath);

    subjectID = config.subjectID;
    sesnum    = config.sessionNum;

    % Support both 'channel' and 'channels' field names
    alignLoss = lossData.allData.(subjectID).session(sesnum).alignment.loss;
    alignWin  = winData.allData.(subjectID).session(sesnum).alignment.win;
    if isfield(alignLoss, 'channel'); lossStruct = alignLoss.channel; else; lossStruct = alignLoss.channels; end
    if isfield(alignWin, 'channel');  winStruct  = alignWin.channel;  else; winStruct  = alignWin.channels;  end

    lossChs = fieldnames(lossStruct);
    winChs  = fieldnames(winStruct);
    sharedChs = intersect(lossChs, winChs);

    mi_loss = zeros(length(sharedChs), 1);
    mi_win  = zeros(length(sharedChs), 1);
    diffs   = zeros(length(sharedChs), 1);

    for i = 1:length(sharedChs)
        ch = sharedChs{i};
        mi_loss(i) = mean(lossStruct.(ch).finalAggMI);
        mi_win(i)  = mean(winStruct.(ch).finalAggMI);
        diffs(i)   = mi_loss(i) - mi_win(i);
    end

    % Paired t-test
    [h, p, ci, stats] = ttest(mi_loss, mi_win);

    % Wilcoxon signed-rank test
    [p_wilcoxon, h_wilcoxon, stats_wilcoxon] = signrank(mi_loss, mi_win);

    outputDir = fullfile(config.outputRoot, sprintf('%s_session%d_stats', subjectID, sesnum));
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Save table
    T = table(sharedChs, mi_loss, mi_win, diffs, ...
        'VariableNames', {'Channel', 'MI_Loss', 'MI_Win', 'Difference'});
    writetable(T, fullfile(outputDir, 'mi_summary.csv'));

    % Save stats report
    fid = fopen(fullfile(outputDir, 'stats_report.txt'), 'w');
    fprintf(fid, 'Subject: %s | Session: %d\n', subjectID, sesnum);
    fprintf(fid, 'Paired t-test (Loss vs Win)\n');
    fprintf(fid, '  t-statistic: %.4f\n', stats.tstat);
    fprintf(fid, '  df: %d\n', stats.df);
    fprintf(fid, '  p-value: %.4e\n', p);
    if p < 0.05
        fprintf(fid, '  Conclusion: Significant difference in MI between conditions (t-test).\n');
    else
        fprintf(fid, '  Conclusion: No significant difference in MI (t-test).\n');
    end
    fprintf(fid, '\nWilcoxon signed-rank test:\n');
    fprintf(fid, '  Signed-rank statistic (W): %d\n', stats_wilcoxon.signedrank);
    fprintf(fid, '  p-value: %.4e\n', p_wilcoxon);
    if p_wilcoxon < 0.05
        fprintf(fid, '  Conclusion: Significant difference in MI between conditions (Wilcoxon).\n');
    else
        fprintf(fid, '  Conclusion: No significant difference in MI (Wilcoxon).\n');
    end
    fclose(fid);

    % Plot 1: Mean MI Bar Plot
    figure;
    bar([mean(mi_loss), mean(mi_win)]);
    set(gca, 'XTickLabel', {'Loss', 'Win'});
    ylabel('Mean Modulation Index');
    title('Mean PAC by Condition');
    saveas(gcf, fullfile(outputDir, 'mean_mi_comparison.png'));
    close;

    % Plot 2: MI Boxplot
    figure;
    boxplot([mi_loss, mi_win], {'Loss', 'Win'});
    ylabel('Modulation Index');
    title('MI Distribution Across Channels');
    saveas(gcf, fullfile(outputDir, 'mi_boxplot.png'));
    close;

    % Plot 3: Histogram with Mean Line
    figure;
    histogram(diffs, 20, 'Normalization', 'pdf');
    hold on;
    xline(0, 'r--', 'Zero');
    xline(mean(diffs), 'b-', 'LineWidth', 2);
    title('Histogram of MI Differences (Loss - Win)');
    xlabel('MI Difference');
    ylabel('Density');
    saveas(gcf, fullfile(outputDir, 'mi_difference_histogram.png'));
    close;

    % Plot 4: Paired Line Plot (Spaghetti Plot)
    figure;
    for i = 1:length(mi_loss)
        plot([1 2], [mi_loss(i), mi_win(i)], '-o', 'Color', [0.6 0.6 0.6]);
        hold on;
    end
    plot([1 2], [mean(mi_loss), mean(mi_win)], '-ok', 'LineWidth', 2);
    xlim([0.8 2.2]);
    xticks([1 2]);
    xticklabels({'Loss', 'Win'});
    ylabel('Modulation Index');
    title('Paired MI per Channel');
    saveas(gcf, fullfile(outputDir, 'paired_line_plot.png'));
    close;

    % Plot 5: Scatter Plot with Identity Line
    figure;
    scatter(mi_win, mi_loss, 30, 'filled');
    hold on;
    plot([0 1], [0 1], 'r--');
    xlabel('MI (Win)');
    ylabel('MI (Loss)');
    title('Channel-wise MI: Loss vs Win');
    axis equal;
    saveas(gcf, fullfile(outputDir, 'scatter_mi_win_vs_loss.png'));
    close;

    % Console summary
    fprintf('Mean MI (loss): %.4f\n', mean(mi_loss));
    fprintf('Mean MI (win):  %.4f\n', mean(mi_win));
    fprintf('Mean difference: %.4f\n', mean(diffs));
    fprintf('Std of difference: %.4f\n', std(diffs));
    fprintf('Paired t-test p-value: %.4e | t = %.4f | df = %d\n', p, stats.tstat, stats.df);
    fprintf('Wilcoxon signed-rank test: p = %.4e | signed-rank = %d\n', p_wilcoxon, stats_wilcoxon.signedrank);
    fprintf('Stats and plots saved to: %s\n', outputDir);
end
