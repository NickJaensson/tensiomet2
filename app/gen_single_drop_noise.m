
close all; clear

% physical parameters
params_phys.sigma = 4;        % surface tension
params_phys.grav = 1.0;       % gravitational acceleration
params_phys.rneedle = 1.0;    % radius of the needle
params_phys.volume0 = 16;     % prescribed volume
params_phys.deltarho = 1.0;   % density difference

% numerical parameters
params_num.N = 40;          % resolution of the discretization for calculation
params_num.Nplot = 80;      % resolution of the discretization for plotting
params_num.Ncheb = 10;      % number of Chebyshev to describe the shape
params_num.alpha = 1;       % relaxation parameter in the Newton-Raphson scheme
params_num.eps = 1e-12;     % convergence critertion: rms(u) < eps
params_num.maxiter = 100;   % maximum number of iteration steps
params.cheb_eps = 1e-2;     % error for describing the shape

% solve the Young-Laplace equation for the given parameters
[vars_sol,vars_num] = solve_forward_young_laplace(params_phys, params_num);

% calculate the volume and the area
volume = pi*vars_num.ws*(vars_sol.r.^2.*sin(vars_sol.psi));
area = pi*2*vars_num.ws*(vars_sol.r);

disp(['volume = ', num2str(volume,15)]);
disp(['area = ', num2str(area,15)]);
disp(['pressure = ', num2str(vars_sol.p0,15)]);

% determine the normal vectors
normals(:,1) = vars_num.Ds*vars_sol.z; % r-component of the normals
normals(:,2) = -vars_num.Ds*vars_sol.r; % z-component of the normals
for i=1:size(normals,2)
    normals(i,:) = normals(i,:)/norm(normals(i,:));
end

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(vars_num.s(1),vars_num.s(end),params_num.Nplot)';
rr = interp1(vars_num.s,vars_sol.r,ss,'pchip');
zz = interp1(vars_num.s,vars_sol.z,ss,'pchip');

nnormals(:,1) = interp1(vars_num.s,normals(:,1),ss,'pchip');
nnormals(:,2) = interp1(vars_num.s,normals(:,2),ss,'pchip');
for i=1:size(nnormals,2)
    nnormals(i,:) = nnormals(i,:)/norm(nnormals(i,:));
end

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rr',zz','b');
plot(rr',zz','b');
set(gca,'DataAspectRatio',[1 1 1])
quiver(rr,zz,nnormals(:,1),nnormals(:,2));

% add noise to the data points
rng(1,"twister");
sigma_noise = 0.001*params_phys.rneedle;
tmp=normrnd(0,sigma_noise,[params_num.Nplot,1]);
for i=1:params_num.Nplot
    rr_noise(i) = rr(i) + tmp(i)*nnormals(i,1);
    zz_noise(i) = zz(i) + tmp(i)*nnormals(i,2);
end
figure; hold on
scatter(rr_noise',zz_noise','b');
plot(rr',zz','r');
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
fr = chebfun(rrb',[ssb(1),ssb(end)],'equi','eps',params.cheb_eps);
fz = chebfun(zzb',[ssb(1),ssb(end)],'equi','eps',params.cheb_eps);
figure; plot(fr); hold on; scatter(ssb,rrb);
xlabel('s','FontSize',24); ylabel('r','FontSize',24);
figure; plot(fz); hold on; scatter(ssb,zzb);
xlabel('s','FontSize',24); ylabel('z','FontSize',24);

vars_num_fit = numerical_grid(params_num,[0,ssb(end)]);

rrfit_nag = gridsample(fr,params_num.N,[0,ssb(end)]);
zzfit_nag = gridsample(fz,params_num.N,[0,ssb(end)]);
psifit_nag = atan2(vars_num_fit.D*zzfit_nag,vars_num_fit.D*rrfit_nag);

% calculate the best fitting Laplace shape
[st,press,rrlaplace,zzlaplace] = ...
    makeIso(zzfit_nag,rrfit_nag,psifit_nag,vars_num_fit.D);

disp(['estimated surface tension = ',num2str(st,12)]);

% current output:
% iter 76: rms(u) = 9.431005e-10
% estimated surface tension = 3.98740839926
