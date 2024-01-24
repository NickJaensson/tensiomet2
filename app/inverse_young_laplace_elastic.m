
gen_single_drop_elastic; 

% close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-9;    % convergence critertion forward: rms(u) < eps
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 100; % maximum number of iteration steps inverse


% the following code can be used to test the CMD implementation
[sigmas, sigmap] = makeCMD(params_phys, vars_sol.psi, vars_sol.r, ...
                           vars_sol.z, vars_num, vars_sol.p0);

plot_surface_stress(vars_num.s, sigmas, sigmap, 2);

[moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure,...
    vars_sol_ref,vars_num_ref,vars_sol,vars_num,params_num);

errorG = abs(moduliS(1)-params_phys.Gmod)/params_phys.Gmod;
errorK = abs(moduliS(2)-params_phys.Kmod)/params_phys.Kmod;

disp(['Error in G = ', num2str(errorG*100,4), ' %']);
disp(['Error in K = ', num2str(errorK*100,4), ' %']);

plot_surface_strain(vars_num.s, lambda_s, lambda_r, 3);
