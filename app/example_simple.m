% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% load the parameter values

parameters_numerical;
parameters_simple;

% solve for the droplet shape (Young-Laplace)

[vars_num, vars_sol] = gen_single_drop(params_phys, params_num, true);

% post processing and plotting

[volume, area] = calculate_volume_area(vars_sol, vars_num, true);

[kappas, kappap] = find_curvature(vars_sol, vars_num);

plot_shape(vars_sol.r, vars_sol.z, 1);
plot_curvature(vars_sol.z, kappas, kappap, 2);