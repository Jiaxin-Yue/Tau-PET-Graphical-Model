function generatePatchInfo(suvr_dir, label_dir, output_dir, subject_id)
%% 
% generate patch information from SUVR and Reeb graph patches
%
% Inputs:
% suvr_dir: the path of SUVR data
% label_dir: the path of Reeb graph
% output_dir: output directory of the patch information
% subject_id


%%

    fid = fopen(suvr_dir,'rb');
    suvr = fread(fid, 'double');
    fclose(fid);
    
    fid = fopen(label_dir,'rb');
    L2 = fread(fid, 'int32');
    fclose(fid);
    L2(L2>100)=0;
    
    Patch={};
    for j=1:max(L2)
        P={}; P{1,1} = j; P{1,2} = sum(L2==j); P{1,3} = {find(L2==j)}; P{1,4} = max(suvr(L2==j));
        P{1,5} = find(suvr==max(suvr(L2==j))); P{1,6} = min(suvr(L2==j));
        Patch(j,:) = P;
    end
    
    if ~exist(strcat(output_dir, subject_id),'dir')
        mkdir(strcat(output_dir, subject_id))
    end
    save((strcat(output_dir,subject_id,'/', subject_id,'_Patch_RH.mat')),'Patch')