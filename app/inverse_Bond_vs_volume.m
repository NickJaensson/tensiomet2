% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% load the parameter values

% numerical parameters
params_num.N = 40;          % grid points for calculation
params_num.eps_fw_simple = 1e-12;  % convergence criterion forward
params_num.maxiter_simple = 100;   % maximum number of iteration steps

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-4;      % error for describing the shape

% parameters for Young-Laplace inverse problem
params_num.eps_inv_yl = 1e-8;    % convergence critertion inverse problem
params_num.sigma_guess = 10;     % guess for interfacial tension value
params_num.p0_guess = 5;         % guess for pressure
params_num.alpha_yl = 0.2;       % relaxation parameter in inverse problem
params_num.maxiter_inv_yl = 5000; % maximum number of iteration steps invers
 
rng(1); % set seed for reproducibility

% physical parameters for the simple droplet problem

params_phys.sigma = 1;      % surface tension
params_phys.grav = 1;       % gravitational acceleration
params_phys.rneedle = 1;    % radius of the needle

Bond_all = 0.5:-0.05:0.10;  % Bond number
Nu_all   = 2:2:22;          % dimensionless volume (Nu)

% parameters for generating the artificial interface points

Nsample = 80;  % number of sample points on interface
               % NOTE: total number of points will be 2*Nsample-1
sigma_noise = 1e-2*params_phys.rneedle; % noise level for sampled points

error_mean_all = zeros(length(Bond_all),length(Nu_all));
Wo_mat = zeros(length(Bond_all),length(Nu_all));

for iii = 1:length(Bond_all)
    
    for jjj = 1:length(Nu_all)

        try            
            Nrand = 200; % nunber of random samples

            params_phys.deltarho = Bond_all(iii);  % density difference
            params_phys.volume0 = Nu_all(jjj);   % prescribed volume
            
            % calculate and display the Worthing number
            params_phys.Wo = params_phys.deltarho*params_phys.grav*params_phys.volume0/...
                (2*pi*params_phys.sigma*params_phys.rneedle);

            Wo_mat(iii,jjj) = params_phys.Wo;

            [vars_num, vars_sol] = gen_single_drop(params_phys, ...
                params_num, false);

            vars_sol.normals = get_normals(vars_sol, vars_num);

            for k = 1:Nrand
    
                [rr_noise, zz_noise] = generate_noisy_shape(vars_sol, vars_num, ...
                    Nsample, sigma_noise);
                
                [vars_sol_fit, vars_num_fit] = ...
                    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);
                
                [st{iii,jjj}(k), press, rrlaplace, zzlaplace] = solve_inverse_young_laplace ( ...
                    vars_sol_fit, params_phys, params_num, vars_num_fit, false);

            end

            % surface tension is 1 in the non-dimensionalization used here
            st_error = abs(st{iii,jjj}-1);
            error_mean_all(iii,jjj) = mean(st_error);

        catch

            error_mean_all(iii,jjj) = NaN;

        end
    end
end

% plot the errors with the NaNs removed
figure
imAlpha = ones(size(error_mean_all));
imAlpha(isnan(error_mean_all))=0;
imagesc(Nu_all, Bond_all, error_mean_all,'AlphaData',imAlpha);
colorbar
set(gca,'color',1.0*[1,1,1]);
set(gca,'ColorScale','log');
set(gca,'YDir','normal');

% figure
% imagesc(Nu_all, Bond_all, Wo_mat); hold on
% colorbar
% contour(Nu_all, Bond_all, Wo_mat,[0.4 0.4],'LineWidth',3,'LineStyle','-');
% % set(gca,'ColorScale','log');
% set(gca,'YDir','normal');
