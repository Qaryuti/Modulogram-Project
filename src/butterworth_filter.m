function filtered_data = butterworth_filter(rawData, Fs, freqRange)
    % BUTTERWORTH_FILTER  Apply a bandpass Butterworth filter.
    %   filtered_data = BUTTERWORTH_FILTER(rawData, Fs, freqRange) applies a
    %   2nd-order bandpass filter to RAWDATA. FREQRANGE specifies [low high]
    %   cutoff frequencies in Hz.

    if numel(freqRange) ~= 2
        error('freqRange must be a 2-element vector [low high].');
    end
    lowFreq  = freqRange(1);
    highFreq = freqRange(2);
    if lowFreq <= 0 || highFreq >= Fs/2
        error('Frequency range must be within (0, Nyquist Frequency).');
    end

    nyquist = Fs/2;
    normalized_range = [lowFreq, highFreq] / nyquist;
    [b, a] = butter(2, normalized_range, 'bandpass');
    filtered_data = filtfilt(b, a, rawData);
end
