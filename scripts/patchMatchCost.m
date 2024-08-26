function patchMatchCost(Patch1,Patch2,dInput, outputCost)
% compute the distance between patches from two subjects
% 

A1 = load(Patch1); A2 = load(Patch2);
load(dInput);

A = zeros(size(A1.Patch,1),size(A2.Patch,1));
for i=1:size(A1.Patch,1)
    for j=1:size(A2.Patch,1)
        A(i,j) = d(A1.Patch{i,5},A2.Patch{j,5})*hausdorff(A1.Patch{i,3}{1,1},A2.Patch{j,3}{1,1},d);
    end
end

save(outputCost,'A');
end