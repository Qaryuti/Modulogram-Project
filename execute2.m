% Example configuration and execution script
config.subjects = {'EMU024','EMU001','EMU025','EMU030'};
config.sessionNum = 1;
config.alignments = {'win', 'loss'};
config.dataDir = '/media/Data/Human_Intracranial_MAD';
config.reference = 'Ground';
config.resultsDir = '/home/alabwam1/Desktop/modulogram_V7/trial-run-3';
config.modulatedRange = [30 100];
config.modWaveRanges.Theta = [4 8];
config.bandwidth = 12;
config.numBins = 12;
config.windowParams = struct('window',[0 1], 'excludeOutOfBounds', true);

modulogram_pipeline(config);
