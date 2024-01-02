
gen_single_drop; 
close all

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-2;   % error for describing the shape
params_num.eps_inv = 1e-9;    % convergence critertion forward: rms(u) < eps
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure
params_num.alpha = 0.5;       % relaxation parameter in inverse problem
params_num.maxiter_inv = 100; % maximum number of iteration steps inverse

Nsample = 40;

[s_plot,r_plot,z_plot] = interpolate_solutions(vars_sol, vars_num, Nsample);

normals_plot = get_normals(vars_sol, vars_num, s_plot);

plot_shape(r_plot, z_plot, 1);
quiver(r_plot,z_plot,normals_plot(:,1),normals_plot(:,2));

% add noise to the data points
rng(1); % set seed
sigma_noise = 0.01*params_phys.rneedle;
tmp=normrnd(0,sigma_noise,[Nsample,1]);
for i=1:Nsample
    rr_noise(i) = r_plot(i) + tmp(i)*normals_plot(i,1);
    zz_noise(i) = z_plot(i) + tmp(i)*normals_plot(i,2);
end
figure; hold on
scatter(rr_noise',zz_noise','b');
plot(r_plot',z_plot','r');
set(gca,'DataAspectRatio',[1 1 1])

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
fr = chebfun(rrb',[ssb(1),ssb(end)],'equi','eps',params_num.eps_cheb);
fz = chebfun(zzb',[ssb(1),ssb(end)],'equi','eps',params_num.eps_cheb);
figure; plot(fr); hold on; scatter(ssb,rrb);
xlabel('s','FontSize',24); ylabel('r','FontSize',24);
figure; plot(fz); hold on; scatter(ssb,zzb);
xlabel('s','FontSize',24); ylabel('z','FontSize',24);

vars_num_fit = numerical_grid(params_num,[0,ssb(end)]);

rr_fit = gridsample(fr,params_num.N,[0,ssb(end)]);
zz_fit = gridsample(fz,params_num.N,[0,ssb(end)]);
psi_fit = atan2(vars_num_fit.D*zz_fit,vars_num_fit.D*rr_fit);

% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = solve_inverse_young_laplace(zz_fit, ...
    rr_fit, psi_fit, params_phys, params_num, vars_num_fit);

disp(['estimated surface tension = ',num2str(st,12)]);

% current output:
% iter 32: rms(u) = 7.486974e-10
% estimated surface tension = 8.10429126183