function patch_channels_field(structFile)
%PATCH_CHANNELS_FIELD Ensure modulogram structs have a 'channels' field.
%   patch_channels_field(FILE) loads FILE (containing allData struct) and
%   copies the .channel field to a new .channels field when missing.
%   The modified structure is saved back to FILE.

if nargin < 1 || ~ischar(structFile)
    error('Specify a modulogram struct .mat file');
end

s = load(structFile);
if ~isfield(s, 'allData')
    error('File %s does not contain allData struct', structFile);
end
subjNames = fieldnames(s.allData);
for si = 1:numel(subjNames)
    subj = subjNames{si};
    sessions = s.allData.(subj).session;
    if ~isstruct(sessions)
        continue;
    end
    for sj = 1:numel(sessions)
        alignNames = fieldnames(sessions(sj).alignment);
        for ai = 1:numel(alignNames)
            al = alignNames{ai};
            alignStruct = sessions(sj).alignment.(al);
            if isfield(alignStruct, 'channel') && ~isfield(alignStruct, 'channels')
                sessions(sj).alignment.(al).channels = alignStruct.channel;
            end
        end
    end
    s.allData.(subj).session = sessions;
end
save(structFile, '-struct', 's');
fprintf('Patched channels field in %s\n', structFile);
end
