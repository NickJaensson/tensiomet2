
example_simple; 
parameters_inverse;

% number of points for the synthetic droplet shape (total amount of points
% will be 2*Nsample-1
Nsample = 80;

% noise level
sigma_noise = 0.01*params_phys.rneedle;
rng(1); % set seed for reproducibility


% generate uniform data points with noise
vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise,zz_noise] = ...
    generate_noisy_shape(vars_sol,vars_num,Nsample,sigma_noise);


% fit the noisy shape with Cheby polynomials
[vars_sol_fit,vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);


% perform inverse Young-Laplace problem
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace ( ...
    vars_sol_fit, params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);


plot_shape(rr_noise, zz_noise, 3);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 3);
plot_shape(rrlaplace, zzlaplace, 3);

legend('noisy data', 'Cheby fit', 'Laplace fit');
