function filtered_data = apply_fir_filter(rawData, Fs, freqRange, fir_order)
%APPLY_FIR_FILTER Apply a zero-phase FIR bandpass filter.
%   filtered_data = APPLY_FIR_FILTER(rawData, Fs, freqRange, fir_order)
%   designs an FIR bandpass filter of order FIR_ORDER using a Hamming
%   window and applies it with zero-phase filtering.

    if numel(freqRange) ~= 2
        error('freqRange must be a 2-element vector [low high].');
    end
    if nargin < 4 || isempty(fir_order)
        fir_order = 1000;
    end

    nyquist = Fs/2;
    lowFreq = max(freqRange(1), 1);
    highFreq = min(freqRange(2), nyquist - eps);
    if lowFreq >= highFreq
        error('Invalid frequency range: [%f %f]', lowFreq, highFreq);
    end

    normalized_range = [lowFreq, highFreq] / nyquist;
    b = fir1(fir_order, normalized_range, 'bandpass', hamming(fir_order + 1));
    filtered_data = filtfilt(b, 1, rawData);
end
