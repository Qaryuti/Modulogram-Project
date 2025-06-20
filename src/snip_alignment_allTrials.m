function [snips, valid_idxs] = snip_alignment_allTrials(data, align_times, Fs, params)
    % SNIP_ALIGNMENT_ALLTRIALS  Return a cell array of trial snippets.
    %   [SNIPS, VALID_IDXS] = SNIP_ALIGNMENT_ALLTRIALS(...) calls SNIP_ALIGNMENT
    %   and converts valid columns to a cell array.

    snipped_matrix = snip_alignment(data, align_times, Fs, params);
    valid_idxs = find(~all(isnan(snipped_matrix),1));
    snips = cell(1, numel(valid_idxs));
    for i = 1:numel(valid_idxs)
        snips{i} = snipped_matrix(:, valid_idxs(i));
    end
end
