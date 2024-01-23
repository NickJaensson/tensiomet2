
gen_single_drop_elastic; 

% close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-9;    % convergence critertion forward: rms(u) < eps
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 0.5;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 1000; % maximum number of iteration steps inverse


% the following code can be used to test the CMD implementation
[sigmas, sigmap] = makeCMD(params_phys, vars_sol.psi, vars_sol.r, ...
                           vars_sol.z, vars_num, vars_sol.p0);

plot_surface_stress(vars_num.s, sigmas, sigmap, 2);


% this step is needed since the makeSFE routine assumes
% a value of C=1
vars_num = numerical_grid(params_num,[0,vars_num.s(end)]);
dummy.C = 1;
vars_num = update_numerical_grid(dummy, vars_num, false);

vars_num_ref = numerical_grid(params_num,[0,vars_num_ref.s(end)]);
dummy.C = 1;
vars_num_ref = update_numerical_grid(dummy, vars_num_ref, false);

global g_strainmeasure g_memptr g_error g_echo glob_w glob_d glob_s 
global glob_r glob_z glob_ts glob_tr

g_strainmeasure = params_phys.strainmeasure;
g_memptr = [1,2];
g_echo = 1;

glob_w = zeros(2,params_num.N);
glob_d = zeros(params_num.N,params_num.N,2);
glob_s = zeros(params_num.N,2);
glob_r = glob_s;
glob_z = glob_s;
glob_ts = glob_s;
glob_tr = glob_s;

% BELOW SHOULD ALL BE CHECKED TO MAKE SURE WE USE THE CORRECT NUMERICAL
% VARIABLES !

glob_w(1,:) = vars_num_ref.ws;
glob_d(:,:,1) = vars_num_ref.Ds;
glob_s(:,1) = vars_num_ref.s;
glob_r(:,1) = vars_sol_ref.r;
glob_z(:,1) = vars_sol_ref.z;
glob_ts(:,1) = params_phys.sigma*ones(vars_num.N,1);
glob_tr(:,1) = params_phys.sigma*ones(vars_num.N,1);

glob_w(2,:) = vars_num.ws;
glob_d(:,:,2) = vars_num.Ds;
glob_s(:,2) = vars_num.s;
glob_r(:,2) = vars_sol.r;
glob_z(:,2) = vars_sol.z;
glob_ts(:,2) = vars_sol.sigmas;
glob_tr(:,2) = vars_sol.sigmap;

[moduliS, lambda_s, lambda_r]  = makeSFE();

errorG = abs(moduliS(1)-params_phys.Gmod)/params_phys.Gmod;
errorK = abs(moduliS(2)-params_phys.Kmod)/params_phys.Kmod;

disp(['Error in G = ', num2str(errorG*100,4), ' %']);
disp(['Error in K = ', num2str(errorK*100,4), ' %']);