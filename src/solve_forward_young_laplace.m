function [vars_sol, vars_num] = solve_forward_young_laplace(params_phys, params_num, shape_guess, vars_num, verbose)
% SOLVE_FORWARD_YOUNG_LAPLACE Solves the Young-Laplace equation for the
% shape of a drop.
%
% This function solves the Young-Laplace equation for a constant surface 
% tension, using a Newton-Raphson iteration.
%
% Usage:
%   [vars_sol, vars_num] = solve_forward_young_laplace(params_phys, params_num)
%
% Inputs:
%   params_phys - A structure containing physical parameters:
%                 deltarho (density difference), grav (gravity), volume0 (initial volume),
%                 sigma (surface tension), rneedle (radius of the needle), 
%                 Wo (Worhtington number).
%   params_num - A structure containing numerical parameters:
%                N (number of points), eps_fw (tolerance for convergence),
%                maxiter (maximum iterations), alpha (relaxation parameter).
%
% Outputs:
%   vars_sol - A structure containing solution variables:
%              r, z (coordinates of the interface), psi (angle of the tangent to the interface),
%              C (stretch parameter), p0 (pressure).
%   vars_num - A structure containing numerical variables used in the solution:
%              Differentiation and integration matrices, Chebyshev points, and modified
%              matrices post-solution (ws, Ds).
%
    
    % interpolate the shape in the Chebyshev points
    r = interp1(shape_guess.s,shape_guess.r,vars_num.s0);
    z = interp1(shape_guess.s,shape_guess.z,vars_num.s0);
    
    psi = atan2(vars_num.D*z,vars_num.D*r);   % intial psi value 
    C = 1;                                % initial stretch parameter
    p0 = 2*params_phys.sigma/params_phys.rneedle;   % initial pressure
    u = ones(3*params_num.N+2,1);             % initial solution vector
    
    % store some variables for the iteration
    iter = 0;
    vars_sol.r = r; vars_sol.z = z; vars_sol.psi = psi;
    vars_sol.C = C; vars_sol.p0 = p0; 
    
    % start the Newton-Raphson iteration
    while rms(u) > params_num.eps_fw_simple
    
        iter = iter + 1;
        
        if iter > params_num.maxiter_simple
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_simple(params_phys,vars_sol,vars_num);
        
        % solve the system of equations
        u = A\b;
        
        % update variables
        vars_sol.r   = vars_sol.r   + u(1:params_num.N);
        vars_sol.z   = vars_sol.z   + u(params_num.N+1:2*params_num.N);
        vars_sol.psi = vars_sol.psi + u(2*params_num.N+1:3*params_num.N); 
        vars_sol.C   = vars_sol.C   + u(3*params_num.N+1);
        vars_sol.p0  = vars_sol.p0  + u(3*params_num.N+2);    
        
        if verbose
            fprintf('iter %d: rms(u) = %d\n',iter,rms(u));
        end

    end

    if any ( vars_sol.r < -1e-6 )
        error('Negative r-coordinates encountered')
    end

    vars_sol.sigmas = params_phys.sigma*ones(vars_num.N,1);
    vars_sol.sigmap = params_phys.sigma*ones(vars_num.N,1);

end