function [AmplP_matrix, MI_vector, binCenters, freq_array] = compute_modulogram_analysis(pac_signals, Fs, num_bins)
    % COMPUTE_MODULOGRAM_ANALYSIS  Compute amplitude-phase matrix and MI values.
    theta_wave  = pac_signals.theta_wave;
    gamma_waves = pac_signals.constituent_gamma_waves;
    freq_array  = pac_signals.center_freqs;

    theta_phase = angle(hilbert(theta_wave));
    binEdges = linspace(-pi, pi, num_bins+1);
    binCenters = 0.5*(binEdges(1:end-1) + binEdges(2:end));

    nFreq = numel(gamma_waves);
    AmplP_matrix = zeros(nFreq, num_bins);
    MI_vector = zeros(nFreq,1);
    for iF = 1:nFreq
        this_gamma = gamma_waves{iF};
        gamma_ampl = abs(hilbert(this_gamma));
        amplP = zeros(1, num_bins);
        for bIdx = 1:num_bins
            mask = theta_phase >= binEdges(bIdx) & theta_phase < binEdges(bIdx+1);
            amplP(bIdx) = mean(gamma_ampl(mask));
        end
        AmplP_matrix(iF,:) = amplP;

        phi = binCenters;
        X = [ones(num_bins,1), cos(phi(:)), sin(phi(:))];
        y = amplP(:);
        params = X \ y;
        fit_vals = X * params;
        rho_sin = sqrt(mean((y - fit_vals).^2));
        MI_vector(iF) = 1 - rho_sin / std(y);
    end
end
