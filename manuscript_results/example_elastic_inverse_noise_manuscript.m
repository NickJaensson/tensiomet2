% calculate the shape for an elastic interface using surface or volume
% compressions / expansions

close all;
clear

% load the parameter values
parameters_numerical_manuscript;
parameters_simple_manuscript;
parameters_elastic_manuscript;
parameters_inverse_manuscript;

% parameters for generating the artificial interface points
Nsample = 80;  % number of sample points on interface
% NOTE: total number of points will be 2*Nsample-1
dimless_noise = 1e-3;
sigma_noise = dimless_noise*params_phys.rneedle; % noise level for sampled points
num_ChebPolys = 8; % number of Chebyshev modes used
Num_Simuls = 10; % number of independent runs to average

% solve for the reference state and the deformed state
[vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
    params_num, true);
[vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
    params_num, vars_num_ref, vars_sol_ref, true);

% calculate the volume and the area of the forward problem
[vars_sol_ref.volume, vars_sol_ref.area] = calculate_volume_area(vars_sol_ref, vars_num_ref, false);
[vars_sol.volume, vars_sol.area] = calculate_volume_area(vars_sol, vars_num, false);

kk = 1;
for i = 1:1000

    if kk >= Num_Simuls + 1
        vars_sol_fit{Num_Simuls+1}.NumRuns = i;
        break
    end
    disp(['Iteration number ', num2str(i), ', successful iteration number ', num2str(kk)]);

    rng(i); % set seed for reproducibility

    % generate uniform data points with noise
    vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
    [rr_noise_ref,zz_noise_ref] = generate_noisy_shape(vars_sol_ref, vars_num_ref, Nsample, sigma_noise);
    vars_sol.normals = get_normals(vars_sol, vars_num);
    [rr_noise,zz_noise] = generate_noisy_shape(vars_sol, vars_num, Nsample, sigma_noise);

    try

        % fit the noisy shape with Cheby polynomials
        [vars_sol_ref_fit{kk},vars_num_ref_fit{kk}] = fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num,num_ChebPolys);
        vars_sol_ref_fit{kk}.p0 = vars_sol_ref.p0;
        vars_sol_ref_fit{kk}.r_noise = rr_noise_ref;
        vars_sol_ref_fit{kk}.z_noise = zz_noise_ref;

        [vars_sol_fit{kk},vars_num_fit{kk}] = fit_shape_with_chebfun(rr_noise,zz_noise,params_num,num_ChebPolys);
        vars_sol_fit{kk}.p0 = vars_sol.p0;
        vars_sol_fit{kk}.r_noise = rr_noise;
        vars_sol_fit{kk}.z_noise = zz_noise;

        if kk < Num_Simuls
            plot_shape(rr_noise_ref*params_phys.rneedle_dimal,zz_noise_ref*params_phys.rneedle_dimal,10);
            plot_shape(vars_sol_ref.r*params_phys.rneedle_dimal,vars_sol_ref.z*params_phys.rneedle_dimal,10);
            plot_shape(vars_sol_ref_fit{kk}.r*params_phys.rneedle_dimal,vars_sol_ref_fit{kk}.z*params_phys.rneedle_dimal,10);
            plot_shape(rr_noise*params_phys.rneedle_dimal,zz_noise*params_phys.rneedle_dimal,10);
            plot_shape(vars_sol.r*params_phys.rneedle_dimal,vars_sol.z*params_phys.rneedle_dimal,10);
            plot_shape(vars_sol_fit{kk}.r*params_phys.rneedle_dimal,vars_sol_fit{kk}.z*params_phys.rneedle_dimal,10);
        end

        % perform CMD to find the surface stresses
        [vars_sol_ref_fit{kk}.sigmas, vars_sol_ref_fit{kk}.sigmap] = makeCMD(params_phys, vars_sol_ref_fit{kk}, vars_num_ref_fit{kk});
        [vars_sol_fit{kk}.sigmas, vars_sol_fit{kk}.sigmap] = makeCMD(params_phys, vars_sol_fit{kk}, vars_num_fit{kk});

        % calculate the volume and the area
        [vars_sol_ref_fit{kk}.volume, vars_sol_ref_fit{kk}.area] = calculate_volume_area(vars_sol_ref_fit{kk}, vars_num_ref_fit{kk}, false);
        [vars_sol_fit{kk}.volume, vars_sol_fit{kk}.area] = calculate_volume_area(vars_sol_fit{kk}, vars_num_fit{kk}, false);


        % provide guesses for K and G
        params_num.K_guess = (vars_sol_fit{kk}.sigmas(1) - vars_sol_ref_fit{kk}.sigmas(1))/log(vars_sol_fit{kk}.area/vars_sol_ref_fit{kk}.area);        % guess for dilational modulus
        params_num.G_guess = params_num.K_guess;          % guess for shear moduls
        params_phys.K_guess = params_num.K_guess;
        params_phys.K_guess2 = (mean([mean(vars_sol_fit{kk}.sigmas), mean(vars_sol_fit{kk}.sigmap)]) - mean([mean(vars_sol_ref_fit{kk}.sigmas), mean(vars_sol_ref_fit{kk}.sigmap)]) )/log(vars_sol_fit{kk}.area/vars_sol_ref_fit{kk}.area);


        % perform SFE to find the moduli and strains
        % NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
        % the the numerical results are used instead of the Cheby fit (giving a
        % best-case scenario)
        try
            [moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure, vars_sol_ref_fit{kk}, vars_num_ref_fit{kk}, vars_sol_fit{kk}, vars_num_fit{kk}, params_num, true);
            if moduliS(2) > 1000
                continue
            end
            vars_sol_fit{kk}.lamp = lambda_r;
            vars_sol_fit{kk}.lams = lambda_s;
            vars_sol_fit{kk}.Gmod = moduliS(1);
            vars_sol_fit{kk}.Kmod = moduliS(2);
            K_dimal = vars_sol_fit{kk}.Kmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
            disp(['Dimensional K: ', num2str(K_dimal), ' mN/m']);
            kk = kk + 1;
        catch
            warning('Make SFE failed');
            moduliS = [0 0];
            lambda_r = zeros(params_num.N,1); lambda_s = zeros(params_num.N,1);
            vars_sol_fit{kk}.lamp = zeros(params_num.N,1);
            vars_sol_fit{kk}.lams = zeros(params_num.N,1);
            vars_sol_fit{kk}.Gmod = 0;
            vars_sol_fit{kk}.Kmod = 0;
        end

    catch

    end
end

% Find averages, for error quantification. Save as last entry of vars_sol_fit cells:
for kk = 1:Num_Simuls
    K_vec(kk) = vars_sol_fit{kk}.Kmod;
    G_vec(kk) = vars_sol_fit{kk}.Kmod;
    sigmas_apex_vec(kk) = vars_sol_fit{kk}.sigmas(1);
    lambdas_apex_vec(kk) = vars_sol_fit{kk}.lams(1);
    lambdap_apex_vec(kk) = vars_sol_fit{kk}.lamp(1);
    volume_vec(kk) = vars_sol_fit{kk}.volume;
    area_vec(kk) = vars_sol_fit{kk}.area;
end
vars_sol_fit{Num_Simuls+1}.Kmod_Avg = mean(K_vec);
vars_sol_fit{Num_Simuls+1}.Kmod_Std = std(K_vec);
vars_sol_fit{Num_Simuls+1}.Gmod_Avg = mean(G_vec);
vars_sol_fit{Num_Simuls+1}.Gmod_Std = std(G_vec);
vars_sol_fit{Num_Simuls+1}.sigmas_apex_Avg = mean(sigmas_apex_vec);
vars_sol_fit{Num_Simuls+1}.sigmas_apex_Std = std(sigmas_apex_vec);
vars_sol_fit{Num_Simuls+1}.lambdas_apex_Avg = mean(lambdas_apex_vec);
vars_sol_fit{Num_Simuls+1}.lambdas_apex_Std = std(lambdas_apex_vec);
vars_sol_fit{Num_Simuls+1}.lambdap_apex_Avg = mean(lambdap_apex_vec);
vars_sol_fit{Num_Simuls+1}.lambdap_apex_Std = std(lambdap_apex_vec);
vars_sol_fit{Num_Simuls+1}.volume_Avg = mean(volume_vec);
vars_sol_fit{Num_Simuls+1}.volume_Std = std(volume_vec);
vars_sol_fit{Num_Simuls+1}.area_Avg = mean(area_vec);
vars_sol_fit{Num_Simuls+1}.area_Std = std(area_vec);

% post processing (re-make dimensional variables)
vars_sol_ref.r_dimal = vars_sol_ref.r*params_phys.rneedle_dimal;
vars_sol_ref.z_dimal = vars_sol_ref.z*params_phys.rneedle_dimal;
vars_sol_ref.p0_dimal = vars_sol_ref.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol_ref.sigmas_dimal = vars_sol_ref.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref.sigmap_dimal = vars_sol_ref.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_ref.volume_dimal = vars_sol_ref.volume*(params_phys.rneedle_dimal)^3;
vars_sol_ref.area_dimal = vars_sol_ref.area*(params_phys.rneedle_dimal)^2;
vars_sol.r_dimal = vars_sol.r*params_phys.rneedle_dimal;
vars_sol.z_dimal = vars_sol.z*params_phys.rneedle_dimal;
vars_sol.p0_dimal = vars_sol.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
vars_sol.sigmas_dimal = vars_sol.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol.sigmap_dimal = vars_sol.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol.volume_dimal = vars_sol.volume*(params_phys.rneedle_dimal)^3;
vars_sol.area_dimal = vars_sol.area*(params_phys.rneedle_dimal)^2;
params_phys.K_guess_dimal = params_phys.K_guess*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
params_phys.K_guess2_dimal = params_phys.K_guess*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.Kmod_Avg_dimal = vars_sol_fit{Num_Simuls+1}.Kmod_Avg*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.Kmod_Std_dimal = vars_sol_fit{Num_Simuls+1}.Kmod_Std*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.Gmod_Avg_dimal = vars_sol_fit{Num_Simuls+1}.Gmod_Avg*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.Gmod_Std_dimal = vars_sol_fit{Num_Simuls+1}.Gmod_Std*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.sigmas_apex_Avg_dimal = vars_sol_fit{Num_Simuls+1}.sigmas_apex_Avg*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.sigmas_apex_Std_dimal = vars_sol_fit{Num_Simuls+1}.sigmas_apex_Std*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.volume_Avg_dimal = vars_sol_fit{Num_Simuls+1}.volume_Avg*(params_phys.rneedle_dimal)^3;
vars_sol_fit{Num_Simuls+1}.volume_Std_dimal = vars_sol_fit{Num_Simuls+1}.volume_Std*(params_phys.rneedle_dimal)^3;
vars_sol_fit{Num_Simuls+1}.area_Avg_dimal = vars_sol_fit{Num_Simuls+1}.area_Avg*(params_phys.rneedle_dimal)^2;
vars_sol_fit{Num_Simuls+1}.area_Std_dimal = vars_sol_fit{Num_Simuls+1}.area_Std*(params_phys.rneedle_dimal)^2;

for kk = 1:Num_Simuls
    vars_sol_ref_fit{kk}.r_dimal = vars_sol_ref_fit{kk}.r*params_phys.rneedle_dimal;
    vars_sol_ref_fit{kk}.z_dimal = vars_sol_ref_fit{kk}.z*params_phys.rneedle_dimal;
    vars_sol_ref_fit{kk}.p0_dimal = vars_sol_ref_fit{kk}.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
    vars_sol_ref_fit{kk}.sigmas_dimal = vars_sol_ref_fit{kk}.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
    vars_sol_ref_fit{kk}.sigmap_dimal = vars_sol_ref_fit{kk}.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
    vars_sol_ref_fit{kk}.volume_dimal = vars_sol_ref_fit{kk}.volume*(params_phys.rneedle_dimal)^3;
    vars_sol_ref_fit{kk}.area_dimal = vars_sol_ref_fit{kk}.area*(params_phys.rneedle_dimal)^2;
    vars_sol_ref_fit{kk}.r_noise_dimal = vars_sol_ref_fit{kk}.r_noise*params_phys.rneedle_dimal;
    vars_sol_ref_fit{kk}.z_noise_dimal = vars_sol_ref_fit{kk}.z_noise*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.r_dimal = vars_sol_fit{kk}.r*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.z_dimal = vars_sol_fit{kk}.z*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.p0_dimal = vars_sol_fit{kk}.p0*params_phys.deltarho_dimal*params_phys.grav_dimal*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.sigmas_dimal = vars_sol_fit{kk}.sigmas*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
    vars_sol_fit{kk}.sigmap_dimal = vars_sol_fit{kk}.sigmap*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
    vars_sol_fit{kk}.volume_dimal = vars_sol_fit{kk}.volume*(params_phys.rneedle_dimal)^3;
    vars_sol_fit{kk}.area_dimal = vars_sol_fit{kk}.area*(params_phys.rneedle_dimal)^2;
    vars_sol_fit{kk}.r_noise_dimal = vars_sol_fit{kk}.r_noise*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.z_noise_dimal = vars_sol_fit{kk}.z_noise*params_phys.rneedle_dimal;
    vars_sol_fit{kk}.Gmod_dimal = vars_sol_fit{kk}.Gmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
    vars_sol_fit{kk}.Kmod_dimal = vars_sol_fit{kk}.Kmod*params_phys.deltarho_dimal*params_phys.grav_dimal*(params_phys.rneedle_dimal)^2;
end
disp(['Dimensional K from Inverse problem: ', num2str(vars_sol_fit{Num_Simuls+1}.Kmod_Avg_dimal), ...
    ' mN/m (STD: ', num2str(vars_sol_fit{Num_Simuls+1}.Kmod_Std_dimal), ' mN/m)']);

%  save data
savefolder = strcat('Results_Noise/Wo=',num2str(params_phys.Wo_paper),'_Ar=',num2str(params_phys.Ar_paper),'_sigma0=',num2str(params_phys.sigma_dimal),...
    '_K=',num2str(params_phys.Kmod_dimal), '_G=',num2str(params_phys.Gmod_dimal),'_StrainM=',params_phys.strainmeasure,'_Strain=', num2str(params_phys.frac),...
    '_Noise=', num2str(dimless_noise), '_NumChebPolys=', num2str(num_ChebPolys), '/');

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
for kk = 1:Num_Simuls
    plot_surface_stress(vars_num_ref.s./vars_num_ref.s(end), vars_sol_ref.sigmas_dimal, vars_sol_ref.sigmap_dimal, 2, '-');
    plot_surface_stress(vars_num_ref_fit{kk}.s./vars_num_ref_fit{kk}.s(end), vars_sol_ref_fit{kk}.sigmas_dimal, vars_sol_ref_fit{kk}.sigmap_dimal, 2, ':');
    title('Ref.: Fwd. (-) vs Inv. (:)');
    if kk == Num_Simuls
        savefig(fullfile(savefolder, 'Surface_Stress_Ref.fig'));
    end

    plot_surface_stress(vars_num.s, vars_sol.sigmas_dimal, vars_sol.sigmap_dimal, 3, '-');
    plot_surface_stress(vars_num_fit{kk}.s, vars_sol_fit{kk}.sigmas_dimal, vars_sol_fit{kk}.sigmap_dimal, 3, ':');
    title('Def.: Fwd. (-) vs Inv. (:)');
    if kk == Num_Simuls
        savefig(fullfile(savefolder, 'Surface_Stress_Deformed.fig'));
    end

    plot_surface_strain(vars_num.s, vars_sol.lams, vars_sol.lamp, 4, '-');
    plot_surface_strain(vars_num_fit{kk}.s, vars_sol_fit{kk}.lams, vars_sol_fit{kk}.lamp, 4, ':');
    title('Def.: Fwd. (-) vs Inv. (:)');
    if kk == Num_Simuls
        savefig(fullfile(savefolder, 'Lambdas_Deformed.fig'));
    end

    plot_surface_deformation(vars_num.s, vars_sol.lams, vars_sol.lamp, 5, '-');
    plot_surface_deformation(vars_num_fit{kk}.s, vars_sol_fit{kk}.lams, vars_sol_fit{kk}.lamp, 5, ':');
    title('Def.: Fwd. (-) vs Inv. (:)');
    if kk == Num_Simuls
        savefig(fullfile(savefolder, 'Deformation_Fractions_Deformed.fig'));
    end

    plot_shape(vars_sol_ref_fit{kk}.r_noise_dimal, vars_sol_ref_fit{kk}.z_noise_dimal, 6);
    plot_shape(vars_sol_ref_fit{kk}.r_dimal, vars_sol_ref_fit{kk}.z_dimal, 6);
    plot_shape(vars_sol_fit{kk}.r_noise_dimal, vars_sol_fit{kk}.z_noise_dimal, 6);
    plot_shape(vars_sol_fit{kk}.r_dimal, vars_sol_fit{kk}.z_dimal, 6);
    if kk == Num_Simuls
        savefig(fullfile(savefolder, 'Drop_Shapes.fig'));
    end

end
