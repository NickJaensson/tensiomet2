function [d, dd, w, s] = numerical_grid(N, domain)
% NUMERICAL_GRID Generates numerical grids and differentiation matrices 
% using Chebyshev points.
%
% This function creates differentiation matrices and integration weights 
% for a given domain using Chebyshev points. It is designed to work with 
% the Chebfun package.
%
% Usage:
%   [d, dd, w, s] = numerical_grid(N, domain)
%
% Inputs:
%   N - The number of points in the grid.
%   domain - The left and right boundaries of the 1D domain, specified as 
%            a vector [left_boundary, right_boundary].
%
% Outputs:
%   d - The first-order differentiation matrix.
%   dd - The second-order differentiation matrix.
%   w - The integration weights.
%   s - The Chebyshev points within the specified domain.
%
% Requirements:
%   This function requires the Chebfun package. Ensure that 'set_paths.m' 
%   in your environment is configured to include the Chebfun installation 
%   path.
%
% Examples:
%   N = 100;
%   domain = [0, 1];
%   [d, dd, w, s] = numerical_grid(N, domain);

%   diffmat, introw and chebpts are defined in the Chebfun package
    d = diffmat(N,1,domain);
    dd = diffmat(N,2,domain);
    w = introw(N,domain);
    s = chebpts(N,domain);

end