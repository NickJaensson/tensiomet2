% creates empty dataset, for runs that failed (i.e., Ar = 1, sigma = 20, frac = 0.8)

close all;
clear

% load the parameter values
parameters_numerical;
parameters_simple;
parameters_elastic;
parameters_inverse;
rng(1); % set seed for reproducibility


% parameters for generating the artificial interface points
Nsample = 80;  % number of sample points on interface
sigma_noise = 0; %1e-4*params_phys.rneedle; % noise level for sampled points


% solve for the reference state (no solution exists for the deformed state)
[vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
    params_num, true);
[vars_sol_ref.volume, vars_sol_ref.area] = calculate_volume_area(vars_sol_ref, vars_num_ref, false);
vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
[rr_noise_ref,zz_noise_ref] = generate_noisy_shape(vars_sol_ref, vars_num_ref, Nsample, sigma_noise);
[vars_sol_ref_fit,vars_num_ref_fit] = fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num);
vars_sol_ref_fit.p0 = vars_sol_ref.p0;
vars_sol_ref_fit.r_noise = rr_noise_ref;
vars_sol_ref_fit.z_noise = zz_noise_ref;
[vars_sol_ref_fit.sigmas, vars_sol_ref_fit.sigmap] = makeCMD(params_phys, vars_sol_ref_fit, vars_num_ref_fit);
[vars_sol_ref_fit.volume, vars_sol_ref_fit.area] = calculate_volume_area(vars_sol_ref_fit, vars_num_ref_fit, false);


% post processing (re-make dimensional variables)
vars_sol_ref.r_dimal = vars_sol_ref.r*params_phys.rneedle_dimal;
vars_sol_ref.z_dimal = vars_sol_ref.z*params_phys.rneedle_dimal;
vars_sol_ref.p0_dimal = vars_sol_ref.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol_ref.sigmas_dimal = vars_sol_ref.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref.sigmap_dimal = vars_sol_ref.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref.volume_dimal = vars_sol_ref.volume*(params_phys.rneedle_dimal)^3;
vars_sol_ref.area_dimal = vars_sol_ref.area*(params_phys.rneedle_dimal)^2;
vars_sol_ref_fit.r_dimal = vars_sol_ref_fit.r*params_phys.rneedle_dimal;
vars_sol_ref_fit.z_dimal = vars_sol_ref_fit.z*params_phys.rneedle_dimal;
vars_sol_ref_fit.p0_dimal = vars_sol_ref_fit.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol_ref_fit.sigmas_dimal = vars_sol_ref_fit.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref_fit.sigmap_dimal = vars_sol_ref_fit.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref_fit.volume_dimal = vars_sol_ref_fit.volume*(params_phys.rneedle_dimal)^3;
vars_sol_ref_fit.area_dimal = vars_sol_ref_fit.area*(params_phys.rneedle_dimal)^2;
vars_sol_ref_fit.r_noise_dimal = vars_sol_ref_fit.r_noise*params_phys.rneedle_dimal;
vars_sol_ref_fit.z_noise_dimal = vars_sol_ref_fit.z_noise*params_phys.rneedle_dimal;

vars_sol.r = zeros(length(vars_sol_ref.r), 1);
vars_sol.z = zeros(length(vars_sol_ref.z), 1);
vars_sol.r_star = zeros(length(vars_sol_ref.r), 1);
vars_sol.z_star = zeros(length(vars_sol_ref.z), 1);
vars_sol.psi = zeros(length(vars_sol_ref.psi), 1);
vars_sol.C = 0;
vars_sol.p0 = 0;
vars_sol.sigmas = zeros(length(vars_sol_ref.sigmas), 1);
vars_sol.sigmap = zeros(length(vars_sol_ref.sigmap), 1);
vars_sol.lams = zeros(length(vars_sol_ref.sigmas), 1);
vars_sol.lamp = zeros(length(vars_sol_ref.sigmap), 1);
vars_sol.volume = 0;
vars_sol.area = 0;
vars_sol.normals = zeros(size(vars_sol_ref.normals, 1), size(vars_sol_ref.normals, 2));
vars_sol_fit.r = zeros(length(vars_sol_ref_fit.r), 1);
vars_sol_fit.z = zeros(length(vars_sol_ref_fit.z), 1);
vars_sol_fit.psi = zeros(length(vars_sol_ref_fit.psi), 1);
vars_sol_fit.p0 = 0;
vars_sol_fit.sigmas = zeros(length(vars_sol_ref_fit.sigmas), 1);
vars_sol_fit.sigmap = zeros(length(vars_sol_ref_fit.sigmap), 1);
vars_sol_fit.lams = zeros(length(vars_sol_ref_fit.sigmas), 1);
vars_sol_fit.lamp = zeros(length(vars_sol_ref_fit.sigmap), 1);
vars_sol_fit.volume = 0;
vars_sol_fit.area = 0;
vars_sol_fit.r_noise = zeros(length(vars_sol_ref_fit.r_noise), 1);
vars_sol_fit.z_noise = zeros(length(vars_sol_ref_fit.z_noise), 1);
vars_sol_fit.Gmod = 0;
vars_sol_fit.Kmod = 0;
params_phys.K_guess = 0;
params_phys.K_guess2 = 0;

vars_sol.r_dimal = zeros(length(vars_sol_ref.r), 1);
vars_sol.z_dimal = zeros(length(vars_sol_ref.z), 1);
vars_sol.p0_dimal = 0;
vars_sol.sigmas_dimal = zeros(length(vars_sol_ref.sigmas), 1);
vars_sol.sigmap_dimal = zeros(length(vars_sol_ref.sigmas), 1);
vars_sol.volume_dimal = 0;
vars_sol.area_dimal = 0;
vars_sol_fit.r_dimal = zeros(length(vars_sol_ref_fit.r), 1);
vars_sol_fit.z_dimal = zeros(length(vars_sol_ref_fit.z), 1);
vars_sol_fit.p0_dimal = 0;
vars_sol_fit.sigmas_dimal = zeros(length(vars_sol_ref_fit.sigmas), 1);
vars_sol_fit.sigmap_dimal = zeros(length(vars_sol_ref_fit.sigmap), 1);
vars_sol_fit.volume_dimal = 0;
vars_sol_fit.area_dimal = 0;
vars_sol_fit.r_noise_dimal = zeros(length(vars_sol_ref_fit.r_noise), 1);
vars_sol_fit.z_noise_dimal = zeros(length(vars_sol_ref_fit.z_noise), 1);
vars_sol_fit.Gmod_dimal = 0;
vars_sol_fit.Kmod_dimal = 0;
params_phys.K_guess_dimal = 0;
params_phys.K_guess2_dimal = 0;

vars_num.D = zeros(params_num.N, params_num.N);
vars_num.DD = zeros(params_num.N, params_num.N);
vars_num.wmat = zeros(params_num.N, params_num.N);
vars_num.w = zeros(1, params_num.N);
vars_num.s0 = zeros(params_num.N, 1);
vars_num.D0 = zeros(params_num.N, params_num.N);
vars_num.w0 = zeros(1, params_num.N);
vars_num.wmat0 = zeros(params_num.N, params_num.N);
vars_num.N = params_num.N;
vars_num.ws = zeros(1, params_num.N);
vars_num.Ds = zeros(params_num.N, params_num.N);
vars_num.s = zeros(params_num.N, 1);
vars_num.wsmat = zeros(params_num.N, params_num.N);
vars_num.C = 0;
vars_num.wsstarmat = zeros(params_num.N, params_num.N);

vars_num_fit.D = zeros(params_num.N, params_num.N);
vars_num_fit.DD = zeros(params_num.N, params_num.N);
vars_num_fit.wmat = zeros(params_num.N, params_num.N);
vars_num_fit.w = zeros(1, params_num.N);
vars_num_fit.s0 = zeros(params_num.N, 1);
vars_num_fit.D0 = zeros(params_num.N, params_num.N);
vars_num_fit.w0 = zeros(1, params_num.N);
vars_num_fit.wmat0 = zeros(params_num.N, params_num.N);
vars_num.N = params_num.N;
vars_num_fit.ws = zeros(1, params_num.N);
vars_num_fit.Ds = zeros(params_num.N, params_num.N);
vars_num_fit.s = zeros(params_num.N, 1);
vars_num_fit.wsmat = zeros(params_num.N, params_num.N);
vars_num.C = 0;

%  save data
savefolder = strcat('Results/Wo=',num2str(params_phys.Wo_paper),'_Ar=',num2str(params_phys.Ar_paper),'_sigma0=',num2str(params_phys.sigma_dimal),...
    '_K=',num2str(params_phys.Kmod_dimal), '_G=',num2str(params_phys.Gmod_dimal),'_StrainM=',params_phys.strainmeasure,'_Strain=', num2str(params_phys.frac), '/');
mkdir(savefolder);
delete(strcat(savefolder, 'params_num.mat'));
save(strcat(savefolder, 'params_num.mat'), 'params_num');
delete(strcat(savefolder, 'params_phys.mat'));
save(strcat(savefolder, 'params_phys.mat'), 'params_phys');
delete(strcat(savefolder, 'vars_num.mat'));
save(strcat(savefolder, 'vars_num.mat'), 'vars_num');
delete(strcat(savefolder, 'vars_num_fit.mat'));
save(strcat(savefolder, 'vars_num_fit.mat'), 'vars_num_fit');
delete(strcat(savefolder, 'vars_num_ref.mat'));
save(strcat(savefolder, 'vars_num_ref.mat'), 'vars_num_ref');
delete(strcat(savefolder, 'vars_num_ref_fit.mat'));
save(strcat(savefolder, 'vars_num_ref_fit.mat'), 'vars_num_ref_fit');
delete(strcat(savefolder, 'vars_sol.mat'));
save(strcat(savefolder, 'vars_sol.mat'), 'vars_sol');
delete(strcat(savefolder, 'vars_sol_fit.mat'));
save(strcat(savefolder, 'vars_sol_fit.mat'), 'vars_sol_fit');
delete(strcat(savefolder, 'vars_sol_ref.mat'));
save(strcat(savefolder, 'vars_sol_ref.mat'), 'vars_sol_ref');
delete(strcat(savefolder, 'vars_sol_ref_fit.mat'));
save(strcat(savefolder, 'vars_sol_ref_fit.mat'), 'vars_sol_ref_fit');



