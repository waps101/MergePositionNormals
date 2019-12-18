function [Vnew] = mergeMeshNormals(V,F,N,lambda )
%MERGEMESHNORMALS Optimise mesh to better match target normals
%   This function takes an initial mesh and target per-vertex normals as 
%   input and optimises the mesh such that its surface normals better 
%   match the target per-vertex normals. It does so by solving a linear
%   system of equations as described in the paper "Efficiently combining 
%   positions and normals for precise 3D geometry", Nehab et al. 2005.
%
% Inputs:
%    V       - nverts x 3 matrix of vertex positions
%    F       - nfaces x 3 matrix of triangle vertex indices
%    N       - nverts x 3 matrix of per-vertex surface normals
%    lambda  - regularisation weight, larger means Vnew is closer to V
%
% Outputs:
%    Vnew    - nverts x 3 refined vertex positions
%
% Implemention created for the following paper which you may like to cite
% in addition to the original Nehab paper:
%
% Ye Yu and William A. P. Smith. Depth estimation meets inverse rendering 
% for single image novel view synthesis. In Proc. CVMP, 2019.

if nargin<4
    lambda=0.1;
end

nverts = size(V,1);

nfaces = size(F,1);

% v1 compared to edge between v2 and v3
r = [1:nfaces        1:nfaces       1:nfaces        1:nfaces        1:nfaces            1:nfaces];
c = [F(:,3)'         F(:,2)'        F(:,3)'+nverts  F(:,2)'+nverts  F(:,3)'+(2*nverts)  F(:,2)'+(2*nverts)];
s = [N(F(:,1),1)'    -N(F(:,1),1)'  N(F(:,1),2)'    -N(F(:,1),2)'   N(F(:,1),3)'        -N(F(:,1),3)'];

% v2 compared to edge between v1 and v3
r = [r nfaces+(1:nfaces)      nfaces+(1:nfaces)       nfaces+(1:nfaces)        nfaces+(1:nfaces)        nfaces+(1:nfaces)            nfaces+(1:nfaces)];
c = [c F(:,3)'       F(:,1)'        F(:,3)'+nverts  F(:,1)'+nverts  F(:,3)'+(2*nverts)  F(:,1)'+(2*nverts)];
s = [s N(F(:,2),1)'  -N(F(:,2),1)'  N(F(:,2),2)'    -N(F(:,2),2)'   N(F(:,2),3)'        -N(F(:,2),3)'];

% v3 compared to edge between v1 and v2
r = [r 2*nfaces+(1:nfaces)      2*nfaces+(1:nfaces)       2*nfaces+(1:nfaces)        2*nfaces+(1:nfaces)        2*nfaces+(1:nfaces)            2*nfaces+(1:nfaces)];
c = [c F(:,2)'       F(:,1)'        F(:,2)'+nverts  F(:,1)'+nverts  F(:,2)'+(2*nverts)  F(:,1)'+(2*nverts)];
s = [s N(F(:,3),1)'  -N(F(:,3),1)'  N(F(:,3),2)'    -N(F(:,3),2)'   N(F(:,3),3)'        -N(F(:,3),3)'];

A = sparse(r,c,s,3*nfaces,(3*nverts));
A = [A; lambda.*speye(3*nverts,3*nverts)];

b = zeros(3*nfaces,1);
b = [b; lambda.*V(:,1); lambda.*V(:,2); lambda.*V(:,3)];

v = A\b;

Vnew(:,1) = v(1:nverts);
Vnew(:,2) = v(nverts+1:2*nverts);
Vnew(:,3) = v(2*nverts+1:3*nverts);

end