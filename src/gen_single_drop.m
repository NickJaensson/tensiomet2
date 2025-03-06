function [vars_num, vars_sol, params_phys] = ...
    gen_single_drop(params_phys, params_num, verbose)
    % GEN_SINGLE_DROP Generates a single drop shape by solving 
    % the Young-Laplace equation.
    %
    % INPUTS:
    %   params_phys - Structure with physical parameters
    %   params_num  - Structure with numerical parameters
    %   verbose     - Boolean flag for displaying progress
    %
    % OUTPUTS:
    %   vars_num    - Updated numerical grid
    %   vars_sol    - Solution variables (shape, pressure, etc.)
    %   params_phys - Updated physical parameters (including area0)

    shape_guess = guess_shape(params_phys, 1000);
    
    vars_num = numerical_grid(params_num, [0,shape_guess.s(end)]);
    
    vars_sol = solve_forward_young_laplace(params_phys, params_num, ...
                                           shape_guess, vars_num, verbose);
    
    vars_num = update_numerical_grid(vars_sol, vars_num, false);

    % store the converged value of area0 in the parameters
    [~, params_phys.area0] = calculate_volume_area(vars_sol, ...
        vars_num, false);

end
