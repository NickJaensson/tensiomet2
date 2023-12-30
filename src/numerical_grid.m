function [vars_num] = numerical_grid(params_num, domain)
% NUMERICAL_GRID Generates numerical grids and differentiation matrices 
% using Chebyshev points.
%
% This function creates differentiation matrices and integration weights 
% for a given domain using Chebyshev points. It is designed to work with 
% the Chebfun package.
%
% Usage:
%   vars_num = numerical_grid(params_num, domain)
%
% Inputs:
%   params_num - A structure containing the number of points in the grid.
%                It should have a field 'N' representing this number.
%   domain - The left and right boundaries of the 1D domain, specified as 
%            a vector [left_boundary, right_boundary].
%
% Output:
%   vars_num - A structure containing the following fields:
%              D: The first-order differentiation matrix.
%              DD: The second-order differentiation matrix.
%              w: The integration weights.
%              s: The Chebyshev points within the specified domain.
%              N: The number of points in the grid (copied from params_num.N).
%
% Requirements:
%   This function requires the Chebfun package. Ensure that 'set_paths.m' 
%   in your environment is configured to include the Chebfun installation 
%   path.

    % diffmat, introw and chebpts are defined in the Chebfun package
    vars_num.D = diffmat(params_num.N,1,domain);
    vars_num.DD = diffmat(params_num.N,2,domain);
    vars_num.w = introw(params_num.N,domain);
    vars_num.s = chebpts(params_num.N,domain);

    vars_num.N = params_num.N; % copy for convenvience

end