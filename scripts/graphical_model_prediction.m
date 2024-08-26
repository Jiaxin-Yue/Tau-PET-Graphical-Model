function subtypes=graphical_model_prediction(subtypes_model,P_training, P, dHausdorff, K, thd, alpha)
% predict the subtype membership for the new data


N = length(subtypes_model);

D = ones(1,N)*inf;
x=1;
for y=1:N
    % compare the temporal orders
    [~,sortedindexY] = sort(cell2mat(P_training{y,1}(:,4)),'descend'); 
    matchIdx1 = findPatchPairs(dHausdorff{x,y}(:,sortedindexY),thd,thd);
    matchIdx1(:,2) = sortedindexY(matchIdx1(:,2));

    [~,sortedindexX] = sort(cell2mat(P(:,4)),'descend'); 
    matchIdx2 = findPatchPairs(dHausdorff{x,y}(sortedindexX,:)',thd,thd);
    matchIdx2(:,2) = sortedindexX(matchIdx2(:,2));
    
      if  sum(ismember(matchIdx1,[sortedindexX(1),sortedindexY(1)],'rows')) || sum(ismember(matchIdx2,[sortedindexY(1),sortedindexX(1)],'rows'))
        % y is taken as the reference
           
        c=0;
        matchIdx1 = findPatchPairs(dHausdorff{x,y}(:,sortedindexY),thd,thd);
        matchIdx1(:,2) = sortedindexY(matchIdx1(:,2));
        pairedSUVRx = cell2mat(P(matchIdx1(:,1),4)); % should been sorted as currentYsequence is sorted
        pairedSUVRy = cell2mat(P_training{y,1}(matchIdx1(:,2),4));
        sortedPatchSize = cell2mat(P(matchIdx1(:,1),2));
        unpairedPatchIdx = setdiff([1:size(P,1)],matchIdx1(:,1));
        Dunpaired1 = sum(cell2mat(P(unpairedPatchIdx,4))./max(cell2mat(P(:,4))).*cell2mat(P(unpairedPatchIdx,2)));
        for i=1:length(pairedSUVRx)-1 % crossing penalty
            for j=i+1:length(pairedSUVRx)
                c = c  + (max(0, pairedSUVRx(j)-pairedSUVRx(i))) * (pairedSUVRy(j) < pairedSUVRy(i))/max(cell2mat(P(:,4))) * sortedPatchSize(i);
            end
        end
        Dcrossing1=c;
        d1 = alpha*Dcrossing1 + Dunpaired1;                
        % x is taken as the reference
        c=0;
        
            
        matchIdx1 = findPatchPairs(dHausdorff{x,y}(sortedindexX,:)',thd,thd);
        matchIdx1(:,2) = sortedindexX(matchIdx1(:,2));
        pairedSUVRy = cell2mat(P_training{y,1}(matchIdx1(:,1),4)); 
        pairedSUVRx = cell2mat(P(matchIdx1(:,2),4)); 
        sortedPatchSize = cell2mat(P_training{y,1}(matchIdx1(:,1),2));
        unpairedPatchIdx = setdiff([1:size(P_training{y,1},1)],matchIdx1(:,1));
        Dunpaired2 = sum(cell2mat(P_training{y,1}(unpairedPatchIdx,4))./max(cell2mat(P_training{y,1}(:,4))).*cell2mat(P_training{y,1}(unpairedPatchIdx,2)));
        for i=1:length(pairedSUVRy)-1 % crossing penalty
            for j=i+1:length(pairedSUVRy)
                c = c  + (max(0, pairedSUVRy(j)-pairedSUVRy(i))) *(pairedSUVRx(j)<pairedSUVRx(i))/max(cell2mat(P_training{y,1}(:,4)))* sortedPatchSize(i);
            end
        end
        Dcrossing2=c;
        
        d2 = alpha*Dcrossing2 + Dunpaired2; 
        D(y) = min(d1,d2);

        
    end
end


[~,maxIDX] = mink(D,K);
subtypes = mode(subtypes_model(maxIDX));
