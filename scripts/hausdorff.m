function dH = hausdorff(A, B, d) 
% A, B: two subsets(patches)
% d: geodesic distances between all vertices on cortex

% A=A(A>=0);
% B=B(B>=0);
dH = max(compute_dist(A, B,d), compute_dist(B, A,d));
end

% Compute distance
function dist = compute_dist(A, B,d) 
m = size(A);
n = size(B);
d_vec = [];
D = [];
% dim= size(A, 2); 
for j = 1:m(1)
    
    for k= 1: n(1)
        
        D(k) = d(A(j),B(k));
      
    end
    
    d_vec(j) = min(D); 
      
end
% keyboard
 dist = max(d_vec);
end