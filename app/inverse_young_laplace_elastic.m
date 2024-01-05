
gen_single_drop_elastic; 
% gen_single_drop; 

close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-3;   % error for describing the shape
params_num.eps_inv = 1e-9;    % convergence critertion forward: rms(u) < eps
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 0.5;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 1000; % maximum number of iteration steps inverse

% number of points for the synthetic droplet shape
Nsample = 80;
Nsample_full = 2*Nsample-1;

% noise level
sigma_noise = 0.01*params_phys.rneedle;

rng(1); % set seed

% interpolate the numerical solutions on a uniform grid.
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
s_plot = linspace(vars_num.s(1),vars_num.s(end),Nsample)';
r_plot = interp1(vars_num.s,vars_sol.r,s_plot,'pchip');
z_plot = interp1(vars_num.s,vars_sol.z,s_plot,'pchip');

normals = get_normals(vars_sol, vars_num);

normals_plot(:,1) = interp1(vars_num.s,normals(:,1),s_plot,'pchip');
normals_plot(:,2) = interp1(vars_num.s,normals(:,2),s_plot,'pchip');
for i=1:size(normals_plot,2)
    normals_plot(i,:) = normals_plot(i,:)/norm(normals_plot(i,:));
end

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

% get continuous s around full shape
ds = sqrt((zz_noise(2:end)-zz_noise(1:end-1)).^2+(rr_noise(2:end)-rr_noise(1:end-1)).^2);
ss_noise = zeros(length(rr_noise),1);
for i=1:length(ds)
    ss_noise(1+i) = sum(ds(1:i));
end

% find smallest stepsize and interpolate on smaller grid
dsmin = min(diff(ss_noise));
ssb = linspace(0,ss_noise(end),ceil(0.5*ss_noise(end)/dsmin));
rrb = interp1(ss_noise,rr_noise,ssb,'spline');
zzb = interp1(ss_noise,zz_noise,ssb,'spline');

% fit the points using Chebyshev polynomials
frr = chebfun(rrb',[ssb(1),ssb(end)],'equi','eps',params_num.eps_cheb);
fzz = chebfun(zzb',[ssb(1),ssb(end)],'equi','eps',params_num.eps_cheb);

coeffs_r = chebcoeffs(frr);
for i = 1:size(coeffs_r,1)
    if rem(i,2) ~= 0 % odd number
        coeffs_r(i) = 0;
    end
end

coeffs_z = chebcoeffs(fzz);
for i = 1:size(coeffs_z,1)
    if rem(i,2) == 0 % even number
        coeffs_z(i) = 0;
    end
end

fr = chebfun(coeffs_r,[ssb(1),ssb(end)],'coeffs');
fz = chebfun(coeffs_z,[ssb(1),ssb(end)],'coeffs');

figure(2); plot(fr); hold on; scatter(ssb,rrb);
xlabel('s','FontSize',24); ylabel('r','FontSize',24);
figure(2); plot(fz); hold on; scatter(ssb,zzb);
xlabel('s','FontSize',24); ylabel('z','FontSize',24);

rr_fit = gridsample(fr,params_num.N,[ssb(end)/2,ssb(end)]);
zz_fit = gridsample(fz,params_num.N,[ssb(end)/2,ssb(end)]);

plot_shape(rr_noise, zz_noise, 3);
plot_shape(rr_fit, zz_fit, 3);

% we use the length of the FITTED shape for the new numerical domain
new_length = integral(sqrt(diff(fr)^2+diff(fz)^2))/2;

% now the mesh is for half of the domain
vars_num_fit = numerical_grid(params_num,[0,new_length]);
dummy.C = 1;
vars_num_fit = update_numerical_grid(dummy, vars_num_fit, false);

psi_fit = atan2(vars_num_fit.Ds*zz_fit,vars_num_fit.Ds*rr_fit);

% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace(zz_fit, ...
    rr_fit, psi_fit, params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);

plot_shape(rr_noise, zz_noise, 4);
plot_shape(rrlaplace, zzlaplace, 4);

% current output:
% iter 195: rms(u) = 9.940579e-10
% estimated surface tension = 3.03343162445

vars_sol_fit.r = rr_fit;
vars_sol_fit.z = zz_fit;
vars_sol_fit.psi = psi_fit;

[kappas,kappap] = find_curvature(vars_sol_fit, vars_num_fit);

plot_curvature(vars_sol.z, kappas, kappap, 5);

[sigmas, sigmap] = makeCMD(params_phys, psi_fit, rr_fit, ...
                           zz_fit, vars_num_fit, vars_sol.p0);

plot_surface_stress(vars_num.s, sigmas, sigmap, 6);

if isfield(vars_sol,'sigmas')
    plot_surface_stress(vars_num.s, vars_sol.sigmas, vars_sol.sigmap, 6);
else
    sigma_vec = params_phys.sigma*ones(vars_num.N,1);
    plot_surface_stress(vars_num.s, sigma_vec, sigma_vec, 6);
end


global g_strainmeasure g_memptr g_error g_echo glob_w glob_d glob_s 
global glob_r glob_z glob_ts glob_tr

g_strainmeasure = 'hencky';
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