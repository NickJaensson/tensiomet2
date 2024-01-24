
% gen_single_drop_elastic; 
gen_single_drop; 

close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-3;    % convergence critertion inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 100; % maximum number of iteration steps inverse

% number of points for the synthetic droplet shape (total amount of points
% will be 2*Nsample-1
Nsample = 80;

% noise level
sigma_noise = 0.01*params_phys.rneedle;

rng(1); % set seed

vars_sol.normals = get_normals(vars_sol, vars_num);

plot_shape(vars_sol.r, vars_sol.z, 1)
quiver(vars_sol.r, vars_sol.z, vars_sol.normals(:,1), vars_sol.normals(:,2));

[rr_noise,zz_noise] = generate_noisy_shape(vars_sol,vars_num,Nsample,sigma_noise);

[vars_sol_fit,vars_num_fit] = fit_shape_with_chebfun(rr_noise,zz_noise,params_num);

plot_shape(rr_noise, zz_noise, 2);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 2);


% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace(vars_sol_fit, ...
    params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);