% load your processed sessStruct
load('analysis/Data/EMU024_session1_sessStruct.mat');

% get channel list
channels = sessStruct.channels;

% collect all anatomical region names
regions = cellfun(@(c) string(c.anatomicalRegion), channels);

% get unique regions
uniqueRegions = unique(regions);

% print
fprintf('Total channels: %d\n', numel(channels));
fprintf('Unique anatomical regions:\n');
for i = 1:numel(uniqueRegions)
    fprintf('  %s\n', uniqueRegions(i));
end
