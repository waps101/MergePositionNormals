# MergePositionNormals

In 2005, Nehab et al. proposed a very efficient method for combining either a depth map or mesh with target surface normals:

    Nehab, D., Rusinkiewicz, S., Davis, J. and Ramamoorthi, R., 2005. Efficiently combining positions and normals for precise 3D geometry. ACM transactions on graphics (TOG), 24(3), pp. 536-543.

The original paper came with a C implementation but this has become more difficult to compile over the years due to required dependencies. In this repository, we provide an efficient re-implementation of the Nehab method in Matlab.

## mergeDepthNormals.m

This merges a perspective depth map and surface normal map.

Inputs:
    Z       - H x W perspective depth map
    NM      - H x W x 3 normal map in cameras coordinates
    cx,cy,f - camera parameters
    lambda  - regularisation weight, larger means Z_merged is closer to Z
    mask    - H x W binary foreground mask (logical)

Outputs:
    Z_merged - H x W refined depth map
    mask     - H x W mask (mask may need to be modified to compute finite differences)

## mergeMeshNormals.m

This merges a mesh and per-vertex surface normals.

# Citation

This implementation was developed for the following paper:

    Ye Yu and William A. P. Smith. Depth estimation meets inverse rendering for single image novel view synthesis. In Proc. CVMP, 2019.

in which we show how to write the depth/normal merging in terms of finite difference matrices and also include the centre of projection. If you use our implementation, you may also like to cite our paper.
