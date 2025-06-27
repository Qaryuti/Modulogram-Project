function render_single_from_struct_file(structFile, subjectID, channelLabel, outDir)
    % Load the .mat file
    s = load(structFile);

    % Navigate to the correct nested struct
    try
        channelField = ['CH' channelLabel];  % e.g., 'CH057'
        alignWin = s.allData.(subjectID).session.alignment.win;
        if isfield(alignWin, 'channel')
            resultStruct = alignWin.channel.(channelField);
        elseif isfield(alignWin, 'channels')
            resultStruct = alignWin.channels.(channelField);
        else
            error('No channel data found');
        end
    catch
        error('Channel %s not found in struct at expected location.', channelField);
    end

    % Extract data for plotting
    AmplP = resultStruct.finalAggAmplP;
    FREQS = resultStruct.centerFreqs;
    PHASES = resultStruct.binCenters;
    region = resultStruct.anatomicalRegion;

    % Generate surface
    [PHASE, FREQ] = meshgrid(PHASES, FREQS);
    fig = figure('Visible','off');
    surf(PHASE, FREQ, AmplP, 'EdgeColor', 'none');
    xlabel('\theta phase (rad)');
    ylabel('\gamma frequency (Hz)');
    zlabel('Normalized amplitude');
    title(sprintf('%s – %s', channelField, region));
    view([45 25]);
    colormap turbo;
    shading interp;
    colorbar;

    % Save output
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    [~, base, ~] = fileparts(structFile);
    outName = fullfile(outDir, sprintf('%s_%s_Surface.png', base, channelField));
    saveas(fig, outName);
    fprintf('✅ Saved surface plot to: %s\n', outName);
end
