function [vars_num] = numerical_grid(params_num, domain)
    % NUMERICAL_GRID Generates numerical grids and differentiation 
    % matrices using Chebyshev points.
    %
    % This function creates differentiation matrices and integration 
    % weights for a given domain using Chebyshev points. 
    % NOTE: this function requires the Chebfun library
    %
    % Inputs:
    %   params_num - Struct with numerical parameters (only parameter used 
    %                here is 'N': number of grid points).
    %   domain     - [left_boundary, right_boundary] of the 1D domain.
    %
    % Output:
    %   vars_num   - Struct with fields:
    %                D   - First-order differentiation matrix.
    %                DD  - Second-order differentiation matrix.
    %                w   - Integration weights.
    %                s   - Chebyshev points in the domain.
    %                N   - Number of grid points.

    % diffmat, introw and chebpts are defined in the Chebfun package
    vars_num.D = diffmat(params_num.N,1,domain);
    vars_num.DD = diffmat(params_num.N,2,domain);
    vars_num.wmat = intmat(params_num.N,1,domain);
    vars_num.w = introw(params_num.N,domain);
    vars_num.s0 = chebpts(params_num.N,domain);

    vars_num.D0 = vars_num.D;
    vars_num.w0 = vars_num.w;
    vars_num.wmat0 = vars_num.wmat;

    vars_num.N = params_num.N; % copy for convenvience

end
