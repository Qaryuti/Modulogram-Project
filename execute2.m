% Example configuration and execution script
config.subjects = {'EMU024','EMU001','EMU025','EMU030','EMU036','EMU037',...
    'EMU038','EMU039','EMU040','EMU041','EMU047','EMU051'};
config.sessionNum = 1;
config.alignments = {'win'};
config.dataDir = '/media/Data/Human_Intracranial_MAD';
config.reference = 'Ground';
config.resultsDir = '/home/alabwam1/Desktop/modulogram_V7/filter-2';
config.modulatedRange = [30 130];
config.modWaveRanges.Theta = [4 8];
config.bandwidth = 7.5;
config.numBins = 32;
config.windowParams = struct('window',[0.7 0.7], 'excludeOutOfBounds', true);

modulogram_pipeline(config);
