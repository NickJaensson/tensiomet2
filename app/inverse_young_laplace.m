
gen_single_drop; 

% close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-4;    % convergence critertion inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 100; % maximum number of iteration steps inverse

% number of points for the synthetic droplet shape (total amount of points
% will be 2*Nsample-1
Nsample = 80;

% noise level
sigma_noise = 0.01*params_phys.rneedle;
rng(1); % set seed for reproducibility


% % run the inverse problem on the numerical problem (best case scenario)
% [st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace ( ...
%     vars_sol, params_phys, params_num, vars_num);
% 
% disp(['estimated surface tension = ',num2str(st,12)]);
% 
% %  run CMD on the numerical problem (best case scenario)
% [sigmas, sigmap] = makeCMD(params_phys, vars_sol, vars_num);
% 
% plot_surface_stress(vars_num.s, sigmas, sigmap, 3);


% generate uniform data points with noise
vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise,zz_noise] = ...
    generate_noisy_shape(vars_sol,vars_num,Nsample,sigma_noise);


% fit the noisy shape with Cheby polynomials
[vars_sol_fit,vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);


% perform inverse Young-Laplace problem
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace ( ...
    vars_sol_fit, params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);


plot_shape(rr_noise, zz_noise, 3);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 3);
plot_shape(rrlaplace, zzlaplace, 3);

legend('noisy data', 'Cheby fit', 'Laplace fit');
