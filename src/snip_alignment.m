function snipped_data = snip_alignment(data, align_times, Fs, params)
    % SNIP_ALIGNMENT  Extract segments around each alignment event.
    %   SNIPPED_DATA = SNIP_ALIGNMENT(data, align_times, Fs, params) extracts
    %   windows of data centered on ALIGN_TIMES (samples). PARAMS.window defines
    %   the [pre post] window in seconds.

    if ~isfield(params, 'window')
        error('params.window is required');
    end
    if ~isfield(params, 'excludeNaNTrials')
        params.excludeNaNTrials = true;
    end

    pre_samples  = round(params.window(1) * Fs);
    post_samples = round(params.window(2) * Fs);
    seg_len = pre_samples + post_samples + 1;

    num_align = numel(align_times);
    snipped_data = NaN(seg_len, num_align);
    for i = 1:num_align
        start_idx = align_times(i) - pre_samples;
        end_idx   = align_times(i) + post_samples;
        if start_idx < 1 || end_idx > length(data)
            if params.excludeNaNTrials
                continue;
            end
        else
            snipped_data(:, i) = data(start_idx:end_idx);
        end
    end

    if params.excludeNaNTrials
        snipped_data(:, all(isnan(snipped_data),1)) = [];
    end
end
