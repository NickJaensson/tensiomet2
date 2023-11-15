
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
[itervars,params] = solve_young_laplace(params);

% calculate the volume and the area
volume = pi*params.ws*(itervars.r.^2.*sin(itervars.psi));
area = pi*2*params.ws*(itervars.r);

disp(['volume = ', num2str(volume,15)]);
disp(['area = ', num2str(area,15)]);
disp(['pressure = ', num2str(itervars.p0,15)]);

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(params.s(1),params.s(end),params.Nplot)';
rr = interp1(params.s,itervars.r,ss,'pchip');
zz = interp1(params.s,itervars.z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rr',zz','b');
plot(rr',zz','b');
set(gca,'DataAspectRatio',[1 1 1])

% determine the curvatures
% NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
% solved here by taking kappap(0) = kappas(0)
kappas = params.Ds*itervars.psi;
kappap = kappas;
kappap(2:end) = sin(itervars.psi(2:end))./itervars.r(2:end);

% plot the curvatures versus the z-coordinate
figure;
plot(itervars.z,kappas,'LineWidth',2); hold on
plot(itervars.z,kappap,'LineWidth',2); hold on
plot(itervars.z,kappap+kappas,'LineWidth',2);
xlabel('z','FontSize',32);
ylabel('\kappa','FontSize',32);
legend('\kappa_s','\kappa_\phi','\kappa_s+\kappa_\phi', ...
    'FontSize',24,'Location','northwest');
xlim([itervars.z(1),0])
ax = gca; ax.FontSize = 24;