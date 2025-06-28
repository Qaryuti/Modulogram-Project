load('analysis/Data/EMU024_session3_sessStruct.mat')
allLabels = [];
for c = 1:numel(sessStruct.channels)
    allLabels = [allLabels; sessStruct.channels{c}.trialLabels];
end
unique(allLabels)
