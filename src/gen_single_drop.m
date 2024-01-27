function [vars_num,vars_sol] = gen_single_drop(params_phys, params_num, verbose)

    shape_guess = guess_shape(params_phys, 1000);
    
    vars_num = numerical_grid(params_num, [0,shape_guess.s(end)]);
    
    vars_sol = solve_forward_young_laplace(params_phys, params_num, ...
                                           shape_guess, vars_num, verbose);
    
    vars_num = update_numerical_grid(vars_sol, vars_num, false);

    % store the converged value of area0 in the parameters
    [~, params_phys.area0] = calculate_volume_area(vars_sol, ...
        vars_num, false);

end