
%gen_single_drop_elastic; 
gen_single_drop; 

close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-3;    % convergence critertion inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 100; % maximum number of iteration steps inverse

% number of points for the synthetic droplet shape
Nsample = 80;
Nsample_full = 2*Nsample-1;

% noise level
sigma_noise = 0.01*params_phys.rneedle;

rng(1); % set seed

normals = get_normals(vars_sol, vars_num);

% interpolate all fields on a uniform grid
[s_plot,r_plot] = interpolate_on_uniform_grid(vars_num,vars_sol.r,Nsample);
[~     ,z_plot] = interpolate_on_uniform_grid(vars_num,vars_sol.z,Nsample);
[~     ,normals_plot(:,1)] = interpolate_on_uniform_grid(vars_num,normals(:,1),Nsample);
[~     ,normals_plot(:,2)] = interpolate_on_uniform_grid(vars_num,normals(:,2),Nsample);

[s_plot_full,r_plot_full,z_plot_full,normals_plot_full] = ...
                               mirror_shape(s_plot,r_plot,z_plot,normals_plot);

plot_shape(r_plot_full, z_plot_full, 1)
quiver(r_plot_full, z_plot_full, normals_plot_full(:,1), normals_plot_full(:,2));

% add noise to the data points
tmp = normrnd(0,sigma_noise,[Nsample_full,1]);
rr_noise = zeros(Nsample_full,1); zz_noise = rr_noise;
for i=1:Nsample_full
    rr_noise(i) = r_plot_full(i) + tmp(i)*normals_plot_full(i,1);
    zz_noise(i) = z_plot_full(i) + tmp(i)*normals_plot_full(i,2);
end

plot_shape(rr_noise, zz_noise, 1);

[vars_sol_fit,vars_num_fit] = fit_shape_with_chebfun(rr_noise,zz_noise,params_num);

plot_shape(rr_noise, zz_noise, 3);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 3);


% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace(vars_sol_fit, ...
    params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);