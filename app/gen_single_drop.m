
close all; clear

% physical parameters
params.sigma = 4;        % surface tension
params.grav = 1.2;       % gravitational acceleration
params.rneedle = 1.4;    % radius of the needle
params.volume0 = 16;     % prescribed volume
params.deltarho = 1.1;   % density difference
params.maxiter = 100;    % maximum number of iteration steps
params.eps = 1e-12;      % convergence critertion: rms(u) < eps

% numerical parameters
params.N = 40;          % resolution of the discretization for calculation
params.Nplot = 80;      % resolution of the discretization for plotting
params.Ncheb = 10;      % number of Chebyshev to describe the shape
params.alpha = 1;       % relaxation parameter in the Newton-Raphson scheme

% solve the Young-Laplace equation for the given parameters
[vars_sol,params] = solve_forward_young_laplace(params);

% calculate the volume and the area
volume = pi*params.ws*(vars_sol.r.^2.*sin(vars_sol.psi));
area = pi*2*params.ws*(vars_sol.r);

disp(['volume = ', num2str(volume,15)]);
disp(['area = ', num2str(area,15)]);
disp(['pressure = ', num2str(vars_sol.p0,15)]);

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(params.s(1),params.s(end),params.Nplot)';
rr = interp1(params.s,vars_sol.r,ss,'pchip');
zz = interp1(params.s,vars_sol.z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rr',zz','b');
plot(rr',zz','b');
set(gca,'DataAspectRatio',[1 1 1])

% determine the curvatures
% NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
% solved here by taking kappap(0) = kappas(0)
kappas = params.Ds*vars_sol.psi;
kappap = kappas;
kappap(2:end) = sin(vars_sol.psi(2:end))./vars_sol.r(2:end);

% plot the curvatures versus the z-coordinate
figure;
plot(vars_sol.z,kappas,'LineWidth',2); hold on
plot(vars_sol.z,kappap,'LineWidth',2); hold on
plot(vars_sol.z,kappap+kappas,'LineWidth',2);
xlabel('z','FontSize',32);
ylabel('\kappa','FontSize',32);
legend('\kappa_s','\kappa_\phi','\kappa_s+\kappa_\phi', ...
    'FontSize',24,'Location','northwest');
xlim([vars_sol.z(1),0])
ax = gca; ax.FontSize = 24;