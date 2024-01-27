
example_elastic;  

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-4;   % error for describing the shape
params_num.eps_inv = 1e-5;    % convergence critertion inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 500; % maximum number of iteration steps inverse

% number of points for the synthetic droplet shape (total amount of points
% will be 2*Nsample-1
Nsample = 80;

% noise level
sigma_noise = 1e-4*params_phys.rneedle;
rng(1); % set seed for reproducibility


% generate uniform data points with noise
vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
[rr_noise_ref,zz_noise_ref] = ...
    generate_noisy_shape(vars_sol_ref,vars_num_ref,Nsample,sigma_noise);

vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise,zz_noise] = ...
    generate_noisy_shape(vars_sol,vars_num,Nsample,sigma_noise);


% fit the noisy shape with Cheby polynomials
[vars_sol_ref_fit,vars_num_ref_fit] = ...
    fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num);
vars_sol_ref_fit.p0 = vars_sol_ref.p0;

[vars_sol_fit,vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);
vars_sol_fit.p0 = vars_sol.p0;


% perform CMD to find the surface stresses
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)
[vars_sol_ref_fit.sigmas, vars_sol_ref_fit.sigmap] = ...
    makeCMD(params_phys, vars_sol_ref_fit, vars_num_ref_fit);

[vars_sol_fit.sigmas, vars_sol_fit.sigmap] = ...
    makeCMD(params_phys, vars_sol_fit, vars_num_fit);

plot_surface_stress(vars_num_ref_fit.s, vars_sol_ref_fit.sigmas, ...
    vars_sol_ref_fit.sigmap, 2);
plot_surface_stress(vars_num_fit.s, vars_sol_fit.sigmas, ...
    vars_sol_fit.sigmap, 2);


% perform SFE to find the moduli and strains
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)
[moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure,...
    vars_sol_ref_fit, vars_num_ref_fit, vars_sol_fit, vars_num_fit, ...
    params_num);

errorG = abs(moduliS(1)-params_phys.Gmod)/params_phys.Gmod;
errorK = abs(moduliS(2)-params_phys.Kmod)/params_phys.Kmod;

disp(['Error in G = ', num2str(errorG*100,4), ' %']);
disp(['Error in K = ', num2str(errorK*100,4), ' %']);

plot_surface_strain(vars_num_fit.s, lambda_s, lambda_r, 3);


plot_shape(rr_noise_ref, zz_noise_ref, 5);
plot_shape(vars_sol_ref_fit.r, vars_sol_ref_fit.z, 5);
plot_shape(rr_noise, zz_noise, 5);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 5);
