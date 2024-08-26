function [subtypes,Connections]=graphical_model_train(N, P, dHausdorff, alpha, unpairedC, thd1, thd2, thdm, K, gamma, iters)
%% 
% graphical model training
%
% Inputs:
% N: total number of training data
% alpha: weights between spatial and temporal penalty
% unpairedC: unpaired cost for the patch matching algorithm
% thd1: threshold for measuring the severity difference between spatial
%       proximity
% thd2: threshold for defining tau positive patch
% thdm: threshold for determing the directionality
% K: parameter for KNN
% gamma: resolution for community detection
% iters: iterations of running clustering; if iters=1, generate the
%        subtypes directly.

%%
if ~exist('P','var') || isempty(P)
    error('Error: The variable P is required.');
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

if ~exist('unpairedC','var') || isempty(unpairedC)
    unpairedC = 500;
end

if ~exist('thd1','var') || isempty(thd1)
    thd1 = 0.5;
end

if ~exist('thd2','var') || isempty(thd2)
    thd2 = 1.3;
end

if ~exist('thdm','var') || isempty(thdm)
    thdm = 1000;
end

if ~exist('K','var') || isempty(K)
    K = floor(N/20);
end

if ~exist('gamma','var') || isempty(gamma)
    gamma = 1;
end

if ~exist('iters','var') || isempty(iters)
    iters = 1;
end

%%
% get the pairs and their crossings

D = ones(N)*inf;

for x=1:N-1
    for y=x+1:length(P)
        [~,sortedindexY] = sort(cell2mat(P{y,1}(:,4)),'descend'); 
        [~,sortedindexX] = sort(cell2mat(P{x,1}(:,4)),'descend'); 
        matchIdx1 = findPatchPairs(dHausdorff{x,y}(:,sortedindexY),unpairedC,unpairedC);
        matchIdx1(:,2) = sortedindexY(matchIdx1(:,2));
        matchIdx2 = findPatchPairs(dHausdorff{y,x}(:,sortedindexX),unpairedC,unpairedC);
        matchIdx2(:,2) = sortedindexX(matchIdx2(:,2));
        
        if sum(ismember(matchIdx1,[sortedindexX(1),sortedindexY(1)],'rows')) || sum(ismember(matchIdx2,[sortedindexY(1),sortedindexX(1)],'rows'))
            
            % subject y is taken as the reference
            c=0;
            matchIdx1 = findPatchPairs(dHausdorff{x,y}(:,sortedindexY),unpairedC,unpairedC);
            matchIdx1(:,2) = sortedindexY(matchIdx1(:,2));
            pairedSUVRx = cell2mat(P{x,1}(matchIdx1(:,1),4)); 
            pairedSUVRy = cell2mat(P{y,1}(matchIdx1(:,2),4)); 
            sortedPatchSize = cell2mat(P{x,1}(matchIdx1(:,1),2));
            unpairedPatchIdx = setdiff([1:size(P{x,1},1)],matchIdx1(:,1));
            Dunpaired = sum(cell2mat(P{x,1}(unpairedPatchIdx,4))/max(cell2mat(P{x,1}(:,4))).*cell2mat(P{x,1}(unpairedPatchIdx,2)));
            for i=1:length(pairedSUVRx)-1 % crossing penalty
                for j=i+1:length(pairedSUVRx)
                     c = c  + (max(0, pairedSUVRx(j)-pairedSUVRx(i))) * ((pairedSUVRy(j)-pairedSUVRy(i))<0)/max(cell2mat(P{x,1}(:,4))) * sortedPatchSize(i);
                  
                end
            end
            Dcrossing=c;
            d1 = alpha*Dcrossing + Dunpaired;
         
            % x is taken as the reference
            c=0;
            matchIdx1 = findPatchPairs(dHausdorff{y,x}(:,sortedindexX),unpairedC,unpairedC);
            matchIdx1(:,2) = sortedindexX(matchIdx1(:,2));
            pairedSUVRy = cell2mat(P{y,1}(matchIdx1(:,1),4)); 
            pairedSUVRx = cell2mat(P{x,1}(matchIdx1(:,2),4));
            sortedPatchSize = cell2mat(P{y,1}(matchIdx1(:,1),2));
            unpairedPatchIdx = setdiff([1:size(P{y,1},1)],matchIdx1(:,1));
            Dunpaired = sum(cell2mat(P{y,1}(unpairedPatchIdx,4))/max(cell2mat(P{y,1}(:,4))).*cell2mat(P{y,1}(unpairedPatchIdx,2)));
            for i=1:length(pairedSUVRy)-1 % crossing penalty
                for j=i+1:length(pairedSUVRy)
                    c = c  + (max(0, pairedSUVRy(j)-pairedSUVRy(i)))*((pairedSUVRx(j)-pairedSUVRx(i))<0)/max(cell2mat(P{y,1}(:,4))) * sortedPatchSize(i);
                end
            end
            Dcrossing=c;
            
            d2 = alpha*Dcrossing + Dunpaired;
            D(x,y) = min(d1,d2); D(y,x) = min(d1,d2);
        end
       
    end
end

%%
%A=1./(D/max(D(~isinf(D)),[],'all')+1);
A = exp(-D/1e3);
A = A - diag(diag(A));
%% consider the directions
diConnections = zeros(N);

for x=1:N-1 
    for y=x+1:N
        [~,sortedindexY] = sort(cell2mat(P{y,1}(:,4)),'descend'); 
        [~,sortedindexX] = sort(cell2mat(P{x,1}(:,4)),'descend');
        % take y as reference
        matchIdx1 = findPatchPairs(dHausdorff{x,y}(:,sortedindexY),unpairedC,unpairedC);
        matchIdx1(:,2) = sortedindexY(matchIdx1(:,2));
        pairedSUVRx = cell2mat(P{x,1}(matchIdx1(:,1),4)); % should been sorted as currentYsequence is sorted
        sortedPatchSizex = cell2mat(P{x,1}(matchIdx1(:,1),2));
        pairedSUVRy = cell2mat(P{y,1}(matchIdx1(:,2),4)); % should been sorted as currentYsequence is sorted
        unpairedPatchIdx = setdiff([1:size(P{x,1},1)],matchIdx1(:,1));
        Munpaired = sum((cell2mat(P{x,1}(unpairedPatchIdx,4))-thd2*ones(length(unpairedPatchIdx),1)).*cell2mat(P{x,1}(unpairedPatchIdx,2)));
        SUVRdiff = pairedSUVRx-pairedSUVRy-thd1; SUVRdiff(SUVRdiff<=0)=0;
        Mcrossing=sum(SUVRdiff.*sortedPatchSizex);
        m1 = Munpaired + Mcrossing;
        % take x as reference
        matchIdx1 = findPatchPairs(dHausdorff{y,x}(:,sortedindexX),unpairedC,unpairedC);
        matchIdx1(:,2) = sortedindexX(matchIdx1(:,2));
        pairedSUVRy = cell2mat(P{y,1}(matchIdx1(:,1),4)); 
        sortedPatchSizey = cell2mat(P{y,1}(matchIdx1(:,1),2));
        pairedSUVRx = cell2mat(P{x,1}(matchIdx1(:,2),4)); 
        sortedPatchSizex = cell2mat(P{x,1}(matchIdx1(:,2),2));
        unpairedPatchIdx = setdiff([1:size(P{y,1},1)],matchIdx1(:,1));
        Munpaired = sum((cell2mat(P{y,1}(unpairedPatchIdx,4))-thd2*ones(length(unpairedPatchIdx),1)).*cell2mat(P{y,1}(unpairedPatchIdx,2)));
        SUVRdiff = pairedSUVRy-pairedSUVRx-thd1; SUVRdiff(SUVRdiff<=0)=0;
        Mcrossing=sum(SUVRdiff.*sortedPatchSizey);
        m2 = Munpaired + Mcrossing;
    
        if m2 - m1 >thdm
            diConnections(x,y) = A(x,y);
        elseif m1 - m2 >thdm
            diConnections(y,x) = A(y,x);
        else
            diConnections(x,y) = A(x,y); diConnections(y,x) = A(y,x);
        end
    end
end
%%
Connections = zeros(length(diConnections));
for i=1:length(diConnections)
    [~,linkIDX] = maxk(diConnections(i,:),K);
    Connections(i,linkIDX) = diConnections(i,linkIDX);
end
G = digraph(Connections);


%%
if iters ==1
    [subtypes,Q]=community_louvain(Connections,gamma);
else
    Mconsensus = zeros(length(Connections),iters);
    for t=1:iters
        [Mconsensus(:,t),Q]=community_louvain(Connections,gamma);
    end
    Cmatrix = ones(length(Connections));
    for i=1:length(Connections)-1
        for j=i+1:length(Connections)
            Cmatrix(i,j)=sum(Mconsensus(i,:)==Mconsensus(j,:))/iters; Cmatrix(j,i) = Cmatrix(i,j);
        end    
    end
    [subtypes,Q]=community_louvain(Cmatrix,gamma);
end
