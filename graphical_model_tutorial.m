%
addpath('./scripts/')

%% generate patch information from SUVR and Reeb graph
suvr_dir = './examples/sub-1/sub-1_SUVR.raw';
label_dir = './examples/sub-1/sub-1_Reeb.raw';
output_dir = './examples/';
subject_id = 'sub-1';

generatePatchInfo(suvr_dir, label_dir, output_dir, subject_id)
%% compute patch distance between two subjects
Patch1 = './examples/sub-1/sub-1_Patch_RH.mat';
Patch2 = './examples/sub-2/sub-2_Patch_RH.mat';
dInput = './utils/dRH.mat';
outputCost = './examples/sub-2/sub-2_pairwiseMatchCost_RH.mat';

patchMatchCost(Patch1,Patch2,dInput, outputCost) % this is the distance between two subjects
%% build the distance matrix and clustering
load('results/trained_results.mat');
alpha = 2;
unpairedC = 500;
thd1 = 0.5; thd2 = 1.3; thdm = 1000; K = 2; 
gamma = 1;
iters = 10;
[subtypes,Connections]=graphical_model_train(N, P, dHausdorff, input_dir, alpha, unpairedC, thd1, thd2, thdm, K, gamma, iters);
%% predict the subtyping membership for new data
load('./results/trained_results.mat','subtypes','P');
P_training = P;
subtypes_model = subtypes;

Patch = load(strcat('./examples/', 'sub-3','/','sub-3','_Patch_RH.mat'),'Patch');
P = Patch.Patch;

dH = load(strcat('./examples/', 'sub-3','/','sub-3','_TotalMatchCost_RH.mat'));
dHausdorff = dH.C;

thd=500;
K = 15;
alpha = 2;
subtypes=graphical_model_prediction(subtypes_model,P_training, P, dHausdorff, K, thd, alpha);
