function test_loop_sessions_channels()
    % TEST_LOOP_SESSIONS_CHANNELS
    % Loops through sessions and channels, dynamically loading Fs and data.
    % Prints debugging information for session and channel iteration.

    % Hardcoded Parameters
    subject_ID = 'EMU024';
    num_sessions = 4; % Example: Total sessions
    data_base_dir = '/media/Data/Human_Intracranial_MAD/1_formatted'; % Correct base directory
    alignment = 'feedback'; % Assuming 'feedback' alignment

    % Loop through each session
    for sesnum = 1:num_sessions
        fprintf('Processing Session %d/%d...\n', sesnum, num_sessions);

        % Load setup and raw data files
        setup_file = fullfile(data_base_dir, subject_ID, sprintf('%s_MAD_SES%d_Setup.mat', subject_ID, sesnum));
        raw_file = fullfile(data_base_dir, subject_ID, sprintf('%s_MAD_SES%d_Raw.mat', subject_ID, sesnum));

        % Check if files exist
        if ~exist(setup_file, 'file')
            fprintf('Setup file not found: %s\n', setup_file);
            continue;
        end
        if ~exist(raw_file, 'file')
            fprintf('Raw file not found: %s\n', raw_file);
            continue;
        end

        % Load setup data
        load(setup_file, 'elec_ind', 'filters', 'trial_times', 'trial_words');

        % Load raw data (Fs)
        load(raw_file, 'Fs');
        if ~exist('Fs', 'var')
            error('Sampling frequency (Fs) is missing in the raw file for session %d.', sesnum);
        end
        fprintf('Session %d: Sampling Frequency (Fs): %d Hz\n', sesnum, Fs);

        % Validate `elec_ind`
        if exist('elec_ind', 'var') && ~isempty(elec_ind)
            fprintf('Loaded setup file. Number of channels: %d\n', numel(elec_ind));
        else
            fprintf('elec_ind missing or empty in setup file for session %d.\n', sesnum);
            continue;
        end

        % Get alignment times
        [align_times, trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);

        % Remove NaN values
        remove_ind = isnan(align_times);
        align_times(remove_ind) = [];
        trial_numbers(remove_ind) = [];

        % Convert alignment times to samples
        align_times = round(align_times * Fs);
        fprintf('Alignment times processed: %d valid times.\n', numel(align_times));

        % Loop through each channel
        for chnum = 1:numel(elec_ind)
            channel_file = fullfile(data_base_dir, subject_ID, 'separate_channel_files', ...
                                    sprintf('%s_MAD_SES%d_ch%03d.mat', subject_ID, sesnum, elec_ind(chnum)));

            if ~exist(channel_file, 'file')
                fprintf('Channel file not found: %s\n', channel_file);
                continue;
            end

            fprintf('Processing Channel %03d...\n', elec_ind(chnum));

            % Load channel data
            load(channel_file, 'data');

            % Debugging: Print basic data information
            fprintf('Channel %03d data loaded. Data size: %d samples.\n', elec_ind(chnum), length(data));

            % Placeholder: Indicate where preprocessing and further steps would go
            fprintf('Placeholder for preprocessing Channel %03d in Session %d.\n', elec_ind(chnum), sesnum);
        end
    end

    fprintf('\nSession and channel loop complete for subject %s.\n', subject_ID);
end
