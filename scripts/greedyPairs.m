function matchIdx1=greedyPairs(C,THD)
% C: cost matrix
% THD
 [minCost,idx] = min(C,[],2);
 matchIdx1 = find(minCost<=THD);
 matchIdx1 = [matchIdx1, idx(matchIdx1)];

end