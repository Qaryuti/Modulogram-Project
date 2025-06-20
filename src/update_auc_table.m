function update_auc_table(subjectFolder, channelID, auc_val)
    % UPDATE_AUC_TABLE  Maintain a per-subject AUC_ranking.txt file.
    rankFile = fullfile(subjectFolder, 'AUC_ranking.txt');
    if isfile(rankFile)
        fid = fopen(rankFile,'r');
        C = textscan(fid,'%s %f','Delimiter','\t');
        fclose(fid);
        chanList = C{1};
        aucList = C{2};
    else
        chanList = {};
        aucList = [];
    end

    idx = find(strcmp(chanList, channelID),1);
    if isempty(idx)
        chanList{end+1} = channelID; %#ok<AGROW>
        aucList(end+1) = auc_val; %#ok<AGROW>
    else
        aucList(idx) = auc_val;
    end

    [aucListSorted, sortIdx] = sort(aucList, 'descend');
    chanListSorted = chanList(sortIdx);

    fid = fopen(rankFile,'w');
    for k = 1:numel(chanListSorted)
        fprintf(fid,'%s\t%.6f\n', chanListSorted{k}, aucListSorted(k));
    end
    fclose(fid);
end
