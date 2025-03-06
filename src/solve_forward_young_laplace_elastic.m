function [vars_sol, vars_num] = ...
    solve_forward_young_laplace_elastic(vars_sol, params_phys, ...
                                        params_num, vars_num, verbose)
    % SOLVE_FORWARD_YOUNG_LAPLACE_ELASTIC Solves the elastic Young-Laplace 
    % equation for a drop shape
    %
    % INPUTS:
    %   vars_sol    - Initial solution variables.
    %   params_phys - Physical parameters
    %   params_num  - Numerical parameters
    %   vars_num    - Numerical variables
    %   verbose     - Flag to print iteration information
    %
    % OUTPUTS:
    %   vars_sol - Updated structure with solution variables
    %   vars_num - Updated structure with numerical variables

    % determine the current target volume/area
    if params_phys.compresstype == 1
        params_phys.volume = params_phys.volume0*params_phys.frac;
    else
        params_phys.area = params_phys.area0*params_phys.frac;
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*params_num.N+2,1);

    % start the Newton-Raphson iteration
    while rms(u) > params_num.eps_fw_elastic
    
        iter = iter + 1;
        
        if iter > params_num.maxiter_elastic
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_elastic(params_phys,vars_sol,vars_num);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        vars_sol.r   = vars_sol.r + u(1:params_num.N);
        vars_sol.z   = vars_sol.z + u(params_num.N+1:2*params_num.N);
        vars_sol.psi = vars_sol.psi + u(2*params_num.N+1:3*params_num.N);    
        vars_sol.sigmas = vars_sol.sigmas + u(3*params_num.N+1:4*params_num.N);
        vars_sol.sigmap = vars_sol.sigmap + u(4*params_num.N+1:5*params_num.N);
        vars_sol.lams = vars_sol.lams + u(5*params_num.N+1:6*params_num.N);
        vars_sol.lamp = vars_sol.lamp + u(6*params_num.N+1:7*params_num.N);
        vars_sol.p0  = vars_sol.p0 + u(end);

        if verbose
            fprintf('iter %d: rms(u) = %d\n',iter,rms(u));
        end

    end
   
end
