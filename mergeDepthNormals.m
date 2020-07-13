function [Z_merged,mask] = mergeDepthNormals(Z,NM,cx,cy,f,lambda,mask)
%MERGEDEPTHNORMALS Optimise depth map to better match target normals
%   This function takes an initial depth map and a target normal map as 
%   input and optimises the depth map such that its surface normals better 
%   match those in the target normal map. It does so by solving a linear
%   system of equations as described in the paper "Efficiently combining 
%   positions and normals for precise 3D geometry", Nehab et al. 2005.
%
% Inputs:
%    Z       - H x W perspective depth map
%    NM      - H x W x 3 normal map in cameras coordinates
%    cx,cy,f - camera parameters
%    lambda  - regularisation weight, larger means Z_merged is closer to Z
%    mask    - H x W binary foreground mask (logical)
%
% Outputs:
%    Z_merged - H x W refined depth map
%    mask     - H x W mask (mask may need to be modified to compute finite
%               differences)
%
% Implemention created for the following paper which you may like to cite
% in addition to the original Nehab paper:
%
% Ye Yu and William A. P. Smith. Depth estimation meets inverse rendering 
% for single image novel view synthesis. In Proc. CVMP, 2019.

[ Dx,Dy,mask,~ ] = gradMatrices( mask,'Backward' );

[x,y]=meshgrid(1:size(mask,2),size(mask,1):-1:1);

npix = sum(mask(:));

X = sparse(1:npix,1:npix,x(mask)-cx,npix,npix);
Y = sparse(1:npix,1:npix,y(mask)-cy,npix,npix);

% Build tangent vector matrices
Tx = [-1/f.*X -1/f.*speye(npix); ...
      -1/f.*Y sparse([],[],[],npix,npix); ...
      speye(npix) sparse([],[],[],npix,npix)] * [Dx; speye(npix)];
Ty = [-1/f.*X sparse([],[],[],npix,npix); ...
      -1/f.*Y -1/f.*speye(npix); ...
      speye(npix) sparse([],[],[],npix,npix)] * [Dy; speye(npix)];

% Build matrix for performing dot product with target surface normals
Nx = NM(:,:,1);
Ny = NM(:,:,2);
Nz = NM(:,:,3);
N = sparse([1:npix 1:npix 1:npix],1:3*npix,[Nx(mask); Ny(mask); Nz(mask)],npix,3*npix);

% Solve linear system
z = [lambda.*speye(npix); N*Tx; N*Ty]\[lambda.*Z(mask); zeros(2*npix,1)];

% Put vector result back into depth map
Z_merged = NaN(size(mask));
Z_merged(mask) = z;

end

