% calculate the shape for an elastic interface using surface or volume
% compressions / expansions

close all;
clear

% load the parameter values
parameters_numerical_manuscript;
parameters_simple_manuscript;
parameters_elastic_manuscript;
parameters_inverse_manuscript;
rng(1); % set seed for reproducibility


% parameters for generating the artificial interface points
Nsample = 80;  % number of sample points on interface
% NOTE: total number of points will be 2*Nsample-1
sigma_noise = 0; %1e-4*params_phys.rneedle; % noise level for sampled points


% solve for the reference state and the deformed state
[vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
    params_num, true);
[vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
    params_num, vars_num_ref, vars_sol_ref, true);


% calculate the volume and the area of the forward problem
[vars_sol_ref.volume, vars_sol_ref.area] = calculate_volume_area(vars_sol_ref, vars_num_ref, false);
[vars_sol.volume, vars_sol.area] = calculate_volume_area(vars_sol, vars_num, false);


% generate uniform data points with noise
vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
[rr_noise_ref,zz_noise_ref] = generate_noisy_shape(vars_sol_ref, vars_num_ref, Nsample, sigma_noise);
vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise,zz_noise] = generate_noisy_shape(vars_sol, vars_num, Nsample, sigma_noise);


% fit the noisy shape with Cheby polynomials
[vars_sol_ref_fit,vars_num_ref_fit] = fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num);
vars_sol_ref_fit.p0 = vars_sol_ref.p0;
vars_sol_ref_fit.r_noise = rr_noise_ref;
vars_sol_ref_fit.z_noise = zz_noise_ref;

[vars_sol_fit,vars_num_fit] = fit_shape_with_chebfun(rr_noise,zz_noise,params_num);
vars_sol_fit.p0 = vars_sol.p0;
vars_sol_fit.r_noise = rr_noise;
vars_sol_fit.z_noise = zz_noise;

% perform CMD to find the surface stresses
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a
% best-case scenario)
% NOTE: in the reference state, we could also fit the YL equations to find
% the stresses. Uncomment the code below to use that approach
% [st, ~, ~, ~] = solve_inverse_young_laplace (vars_sol_ref_fit, params_phys, params_num, vars_num_ref_fit);
% vars_sol_ref_fit.sigmas = st*ones(vars_num.N,1);
% vars_sol_ref_fit.sigmap = st*ones(vars_num.N,1);
[vars_sol_ref_fit.sigmas, vars_sol_ref_fit.sigmap] = makeCMD(params_phys, vars_sol_ref_fit, vars_num_ref_fit);
[vars_sol_fit.sigmas, vars_sol_fit.sigmap] = makeCMD(params_phys, vars_sol_fit, vars_num_fit);


% calculate the volume and the area
[vars_sol_ref_fit.volume, vars_sol_ref_fit.area] = calculate_volume_area(vars_sol_ref_fit, vars_num_ref_fit, false);
[vars_sol_fit.volume, vars_sol_fit.area] = calculate_volume_area(vars_sol_fit, vars_num_fit, false);

% provide guesses for K and G
params_num.K_guess = (vars_sol_fit.sigmas(1) - vars_sol_ref_fit.sigmas(1))/log(vars_sol_fit.area/vars_sol_ref_fit.area);        % guess for dilational modulus
params_num.G_guess = params_num.K_guess;          % guess for shear moduls
params_phys.K_guess = params_num.K_guess;
params_phys.K_guess2 = (mean([mean(vars_sol_fit.sigmas), mean(vars_sol_fit.sigmap)]) - mean([mean(vars_sol_ref_fit.sigmas), mean(vars_sol_ref_fit.sigmap)]) )/log(vars_sol_fit.area/vars_sol_ref_fit.area);


% perform SFE to find the moduli and strains
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a
% best-case scenario)
try
    [moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure, vars_sol_ref_fit, vars_num_ref_fit, vars_sol_fit, vars_num_fit, params_num, true);
    vars_sol_fit.lamp = lambda_r;
    vars_sol_fit.lams = lambda_s;
    vars_sol_fit.Gmod = moduliS(1);
    vars_sol_fit.Kmod = moduliS(2);
    disp(['Dimensional K: ', num2str(vars_sol_fit.Kmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2), ' mN/m']);
catch
    warning('Make SFE failed');
    moduliS = [0 0];
    lambda_r = zeros(params_num.N,1); lambda_s = zeros(params_num.N,1);
    vars_sol_fit.lamp = zeros(params_num.N,1);
    vars_sol_fit.lams = zeros(params_num.N,1);
    vars_sol_fit.Gmod = 0;
    vars_sol_fit.Kmod = 0;
end


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
vars_sol.r_dimal = vars_sol.r*params_phys.rneedle_dimal;
vars_sol.z_dimal = vars_sol.z*params_phys.rneedle_dimal;
vars_sol.p0_dimal = vars_sol.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol.sigmas_dimal = vars_sol.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol.sigmap_dimal = vars_sol.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol.volume_dimal = vars_sol.volume*(params_phys.rneedle_dimal)^3;
vars_sol.area_dimal = vars_sol.area*(params_phys.rneedle_dimal)^2;
vars_sol_fit.r_dimal = vars_sol_fit.r*params_phys.rneedle_dimal;
vars_sol_fit.z_dimal = vars_sol_fit.z*params_phys.rneedle_dimal;
vars_sol_fit.p0_dimal = vars_sol_fit.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol_fit.sigmas_dimal = vars_sol_fit.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit.sigmap_dimal = vars_sol_fit.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit.volume_dimal = vars_sol_fit.volume*(params_phys.rneedle_dimal)^3;
vars_sol_fit.area_dimal = vars_sol_fit.area*(params_phys.rneedle_dimal)^2;
vars_sol_fit.r_noise_dimal = vars_sol_fit.r_noise*params_phys.rneedle_dimal;
vars_sol_fit.z_noise_dimal = vars_sol_fit.z_noise*params_phys.rneedle_dimal;
vars_sol_fit.Gmod_dimal = vars_sol_fit.Gmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit.Kmod_dimal = vars_sol_fit.Kmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
params_phys.K_guess_dimal = params_phys.K_guess*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
params_phys.K_guess2_dimal = params_phys.K_guess*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;


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


% plotting
errorG = abs(moduliS(1)-params_phys.Gmod)/params_phys.Gmod;
errorK = abs(moduliS(2)-params_phys.Kmod)/params_phys.Kmod;

disp(['Error in G = ', num2str(errorG*100,4), ' %']);
disp(['Error in K = ', num2str(errorK*100,4), ' %']);

plot_surface_stress(vars_num_ref.s./vars_num_ref.s(end), vars_sol_ref.sigmas_dimal, vars_sol_ref.sigmap_dimal, 2, '-');
plot_surface_stress(vars_num_ref_fit.s./vars_num_ref_fit.s(end), vars_sol_ref_fit.sigmas_dimal, vars_sol_ref_fit.sigmap_dimal, 2, ':');
title('Ref.: Fwd. (-) vs Inv. (:)');
savefig(fullfile(savefolder, 'Surface_Stress_Ref.fig'));

plot_surface_stress(vars_num.s, vars_sol.sigmas_dimal, vars_sol.sigmap_dimal, 3, '-');
plot_surface_stress(vars_num_fit.s, vars_sol_fit.sigmas_dimal, vars_sol_fit.sigmap_dimal, 3, ':');
title('Def.: Fwd. (-) vs Inv. (:)');
savefig(fullfile(savefolder, 'Surface_Stress_Deformed.fig'));

plot_surface_strain(vars_num.s, vars_sol.lams, vars_sol.lamp, 4, '-');
plot_surface_strain(vars_num_fit.s, lambda_s, lambda_r, 4, ':');
title('Def.: Fwd. (-) vs Inv. (:)');
savefig(fullfile(savefolder, 'Lambdas_Deformed.fig'));

plot_surface_deformation(vars_num.s, vars_sol.lams, vars_sol.lamp, 5, '-');
plot_surface_deformation(vars_num_fit.s, lambda_s, lambda_r, 5, ':');
title('Def.: Fwd. (-) vs Inv. (:)');
savefig(fullfile(savefolder, 'Deformation_Fractions_Deformed.fig'));

plot_shape(vars_sol_ref_fit.r_noise_dimal, vars_sol_ref_fit.z_noise_dimal, 6);
plot_shape(vars_sol_ref_fit.r_dimal, vars_sol_ref_fit.z_dimal, 6);
plot_shape(vars_sol_fit.r_noise_dimal, vars_sol_fit.z_noise_dimal, 6);
plot_shape(vars_sol_fit.r_dimal, vars_sol_fit.z_dimal, 6);
savefig(fullfile(savefolder, 'Drop_Shapes.fig'));


