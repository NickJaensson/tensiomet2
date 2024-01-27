
function [vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
    params_num, vars_num_ref,vars_sol_ref)

    vars_num = vars_num_ref;
    vars_sol = vars_sol_ref;
    
    % NOTE: at this stage the initial guesses for r, z, psi and p0 are taken
    % from the solution without elasticity. 
    
    % initialize the surface strains and stresses
    % NOTE: in the solution procedure, the dlams and dlamp (Newton-Raphson 
    % update variables) are required to be equal at s=0. For the correct 
    % solution, the initial guess must have equal lams and lamp at s=0
    vars_sol.lamp = ones(params_num.N,1); vars_sol.lams = vars_sol.lamp;
    vars_sol.sigmas = params_phys.sigma*ones(params_num.N,1); 
    vars_sol.sigmap = vars_sol.sigmas;
    
    vars_sol.r_star = vars_sol.r; vars_sol.z_star = vars_sol.z;
        
    [vars_sol,vars_num] = ...
        solve_forward_young_laplace_elastic(vars_sol, params_phys, ...
                                            params_num, vars_num);

    vars_num = update_numerical_grid(vars_sol, vars_num, true);

end