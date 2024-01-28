% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% load the parameter values

parameters_numerical;
parameters_simple;
parameters_inverse;
rng(1); % set seed for reproducibility

% parameters for generating the artificial interface points

Nsample = 80;  % number of sample points on interface
               % NOTE: total number of points will be 2*Nsample-1
sigma_noise = 1e-2*params_phys.rneedle; % noise level for sampled points

% solve for the droplet shape (Young-Laplace)

[vars_num, vars_sol, params_phys] = gen_single_drop(params_phys, params_num, true);

% generate uniform data points with noise

vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise, zz_noise] = generate_noisy_shape(vars_sol, vars_num, ...
    Nsample, sigma_noise);

% fit the noisy shape with Cheby polynomials

[vars_sol_fit, vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);

% perform inverse Young-Laplace problem
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)

[st, press, rrlaplace, zzlaplace] = solve_inverse_young_laplace ( ...
    vars_sol_fit, params_phys, params_num, vars_num_fit, true);

% post processing and plotting

disp(['estimated surface tension = ',num2str(st,12)]);

plot_shape(rr_noise, zz_noise, 3);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 3);
plot_shape(rrlaplace, zzlaplace, 3);

legend('noisy data', 'Cheby fit', 'Laplace fit');