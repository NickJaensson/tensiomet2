function [s_full, r_full, z_full,normals_full] = ...
    mirror_shape(s, r, z, normals)
    % MIRROR_SHAPE Mirrors a shape symmetrically about the z-axis.
    %
    % Given the shape coordinates (s, r, z) and normal vectors, this 
    % function mirrors the shape across the z-axis to create a full 
    % profile.
    %
    % INPUTS:
    %   s        - Arclength coordinates.
    %   r        - Radial coordinates.
    %   z        - Axial coordinates.
    %   normals  - Normal vectors [nx, nz].
    %
    % OUTPUTS:
    %   s_full       - Full mirrored arclength coordinates.
    %   r_full       - Full mirrored radial coordinates.
    %   z_full       - Full mirrored axial coordinates.
    %   normals_full - Full mirrored normal vectors.

    tmp = s-s(end);
    s_full = [tmp(1:end-1)',s']';
    
    r_full = [-flipud(r)',r(2:end)']';
    
    z_full = [flipud(z)',z(2:end)']';
    
    normals_full(:,1) = [-flipud(normals(:,1))',normals(2:end,1)']';
    normals_full(:,2) = [flipud(normals(:,2))',normals(2:end,2)']';

end
