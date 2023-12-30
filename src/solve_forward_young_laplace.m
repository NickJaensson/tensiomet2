function [vars_sol, vars_num] = solve_forward_young_laplace(params_phys, params_num)
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
%                 sigma (surface tension), rneedle (radius of the needle).
%   params_num - A structure containing numerical parameters:
%                N (number of points), eps (tolerance for convergence),
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

    % calculate the Worthinton number
    params_phys.Wo = params_phys.deltarho*params_phys.grav*params_phys.volume0/...
                                        (2*pi*params_phys.sigma*params_phys.rneedle);

    % find an initial guess of the shape
    [r_guess, z_guess, s_guess] = guess_shape(params_phys,1000);
    
    smax = s_guess(end); % the total length of the 1D domain
    
    % get the differentation/integration matrices and the grid
    vars_num = numerical_grid(params_num,[0,smax]);
    
    % interpolate the shape in the Chebyshev points
    r = interp1(s_guess,r_guess,vars_num.s);
    z = interp1(s_guess,z_guess,vars_num.s);
    
    psi = atan2(vars_num.D*z,vars_num.D*r);   % intial psi value 
    C = 1;                                % initial stretch parameter
    p0 = 2*params_phys.sigma/params_phys.rneedle;   % initial pressure
    u = ones(3*params_num.N+2,1);             % initial solution vector
    
    % store some variables for the iteration
    iter = 0;
    vars_sol.r = r; vars_sol.z = z; vars_sol.psi = psi;
    vars_sol.C = C; vars_sol.p0 = p0; 
    
    % start the Newton-Raphson iteration
    while rms(u) > params_num.eps
    
        iter = iter + 1;
        
        if iter > params_num.maxiter
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_simple(params_phys,vars_sol,vars_num);
        
        % solve the system of equations
        u = A\b;
        
        % update variables
        vars_sol.r   = vars_sol.r   + params_num.alpha*u(1:params_num.N);
        vars_sol.z   = vars_sol.z   + params_num.alpha*u(params_num.N+1:2*params_num.N);
        vars_sol.psi = vars_sol.psi + params_num.alpha*u(2*params_num.N+1:3*params_num.N); 
        vars_sol.C   = vars_sol.C   + params_num.alpha*u(3*params_num.N+1);
        vars_sol.p0  = vars_sol.p0  + params_num.alpha*u(3*params_num.N+2);    
        
        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));
    
    end

    % the integration and differentation matrices in the solution state
    vars_num.ws = vars_num.w/vars_sol.C; 
    vars_num.Ds = vars_sol.C*vars_num.D; 
    
end