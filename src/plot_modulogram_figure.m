function plot_modulogram_figure(AmplP_matrix, ~, binCenters, freq_array, outFile, anatomicalRegion)
    % PLOT_MODULOGRAM_FIGURE  Draw modulogram and MI curve and save to file.
    fig = figure('Visible','off','Position',[100,100,800,600]);
    subplot(2,1,1);
    imagesc('XData', binCenters, 'YData', freq_array, 'CData', AmplP_matrix);
    axis xy; colorbar;
    xlabel('\theta phase (rad)');
    ylabel('\gamma centre frequency (Hz)');
    title(sprintf('Modulogram \x2013 %s', anatomicalRegion));

    [nFreq, nBins] = size(AmplP_matrix);
    MI_new = zeros(nFreq,1);
    phi = binCenters(:);
    X = [ones(nBins,1) cos(phi) sin(phi)];
    for iF = 1:nFreq
        y = AmplP_matrix(iF,:).';
        params = X \ y;
        fit_vals = X * params;
        rho_sin = sqrt(mean((y - fit_vals).^2));
        MI_new(iF) = 1 - rho_sin / std(y);
    end

    subplot(2,1,2);
    plot(freq_array, MI_new, '-o', 'LineWidth', 2);
    xlabel('\gamma centre frequency (Hz)');
    ylabel('Modulation Index');
    grid on;
    title('Modulation Index across \gamma sub-bands');

    print(fig, '-dpng', outFile);
    close(fig);

    auc_val = trapz(freq_array, MI_new);
    sessionFolder = fileparts(outFile);
    subjectFolder = fileparts(sessionFolder);
    [~, pngName, ~] = fileparts(outFile);
    channelID = strtok(pngName, '_');
    update_auc_table(subjectFolder, channelID, auc_val);
end
