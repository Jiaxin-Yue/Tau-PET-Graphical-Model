function matchIdx1=findPatchPairs(C,Cunpair, THD)
% C: cost matrix
% Cunpair: unpair cost
matchIdx1 = matchpairs(C,Cunpair); % firstly apply a linear assignment to find the best matched result

% check if unpaired patch also meet the threshold setting
flag1=1;

while flag1
    unpairedPatchIdx = setdiff([1:size(C,1)],matchIdx1(:,1));
    if ~isempty(unpairedPatchIdx)
            matchIdx2 = matchpairs(C(unpairedPatchIdx,:),THD);
            matchIdx2(:,1) = unpairedPatchIdx(matchIdx2(:,1));
            matchIdx1 = [matchIdx1; matchIdx2];
            if isempty(matchIdx2)
                flag1=0;
            end
    else
        flag1=0;
    end
end

end