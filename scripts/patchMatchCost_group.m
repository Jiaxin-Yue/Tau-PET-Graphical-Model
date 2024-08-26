function patchMatchCost_group(thisSubj,subject_list,dInput1, dInput2, fileDir)
load(subject_list)
Patch1 = strcat(fileDir, thisSubj,'/',thisSubj,'_Patch_RH.mat');
if exist(Patch1,'file')
    A1 = load(Patch1); 
    load(dInput1);
    C = cell(1,length(subject_list));
    outputCost = strcat(fileDir, thisSubj,'/',thisSubj,'_TotalMatchCost_RH.mat');
    if ~exist(outputCost,'file')
    for x = 1:length(subject_list)
        subject_id = subject_list(x);
        disp(subject_id)
        %outputCost = strcat(fileDir, thisSubj,'/',subject_id,'_MatchCost_RH.mat');
        %if ~exist(outputCost,'file')
            Patch2 = strcat(fileDir, subject_id,'/',subject_id,'_Patch_RH.mat');
            A2 = load(Patch2);
            A = zeros(size(A1.Patch,1),size(A2.Patch,1));
            for i=1:size(A1.Patch,1)
                for j=1:size(A2.Patch,1)
                    A(i,j) = d(A1.Patch{i,5},A2.Patch{j,5})*hausdorff(A1.Patch{i,3}{1,1},A2.Patch{j,3}{1,1},d);
                end
            end
        C{1,x} = A;
            %save(outputCost,'A');
        %end
    
    end
    save(outputCost,'C');
    end
end

Patch1 = strcat(fileDir, thisSubj,'/',thisSubj,'_Patch_LH.mat');
if exist(Patch1,'file')
    A1 = load(Patch1);
    load(dInput2);
    C = cell(1,length(subject_list));
    outputCost = strcat(fileDir, thisSubj,'/',thisSubj,'_TotalMatchCost_LH.mat');
    if ~exist(outputCost,'file')
    for x = 1:length(subject_list)
        subject_id = subject_list(x);
        disp(subject_id)
        %outputCost = strcat(fileDir, thisSubj,'/',subject_id,'_MatchCost_LH.mat');
        %if ~exist(outputCost,'file')
            Patch2 = strcat(fileDir, subject_id,'/',subject_id,'_Patch_LH.mat');
            A2 = load(Patch2);
            A = zeros(size(A1.Patch,1),size(A2.Patch,1));
            for i=1:size(A1.Patch,1)
                for j=1:size(A2.Patch,1)
                    A(i,j) = d(A1.Patch{i,5},A2.Patch{j,5})*hausdorff(A1.Patch{i,3}{1,1},A2.Patch{j,3}{1,1},d);
                end
            end
        C{1,x}=A;
            %save(outputCost,'A');
        %end
    end

    save(outputCost,'C');
    end
end
end