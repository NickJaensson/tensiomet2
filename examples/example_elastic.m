% calculate the shape for an elastic interface using surface or volume
% compressions / expansions

close all; clear

% load the parameter values

parameters_numerical;
parameters_simple;
parameters_elastic;

% solve for the reference state and the deformed state

[vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
    params_num, true);

[vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
    params_num, vars_num_ref, vars_sol_ref, true);

% post processing and plotting

[volume, area] = calculate_volume_area(vars_sol, vars_num, true);

plot_shape(vars_sol.r, vars_sol.z, 1);

plot_surface_stress(vars_num.s, vars_sol.sigmas, vars_sol.sigmap, 2);

plot_surface_strain(vars_num.s, vars_sol.lams, vars_sol.lamp, 3);

[kappas, kappap] = find_curvature(vars_sol, vars_num);

plot_curvature(vars_sol.z, kappas, kappap, 4);
