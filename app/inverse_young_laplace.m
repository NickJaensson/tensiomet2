
gen_single_drop; 

close all

% numerical parameters for inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.eps_inv = 1e-9;    % convergence critertion forward: rms(u) < eps
params_num.maxiter_inv = 1000; % maximum number of iteration steps inverse
params_num.alpha = 0.5;       % relaxation parameter in inverse problem


% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace(vars_sol.z, ...
    vars_sol.r, vars_sol.psi, params_phys, params_num, vars_num);

disp(['estimated surface tension = ',num2str(st,12)]);


% the following code can be used to test the CMD implementation
[sigmas, sigmap] = makeCMD(params_phys, vars_sol.psi, vars_sol.r, ...
                           vars_sol.z, vars_num, vars_sol.p0);

plot_surface_stress(vars_num.s, sigmas, sigmap, 3);