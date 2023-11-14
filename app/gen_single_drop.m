
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
volume = pi*params.w*(itervars.r.^2.*sin(itervars.psi))/itervars.C;
area = pi*2*params.w*(itervars.r)/itervars.C;

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