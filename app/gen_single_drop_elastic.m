% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% physical parameters
params.sigma = 4;        % surface tension
params.grav = 1.2;       % gravitational acceleration
params.rneedle = 1.4;    % radius of the needle
params.volume0 = 16;     % prescribed volume
params.deltarho = 1.1;   % density difference
params.maxiter = 100;    % maximum number of iteration steps
params.eps = 1e-12;      % convergence critertion: rms(u) < eps

% physical parameters for the elastic problem
params.Kmod = 3;          % elastic dilational modulus
params.Gmod = 2;          % elastic shear modulus
params.compresstype = 1;  % 1: compress the volume other: compress the area
params.fracm = [0.8];     % compute elastic stresses for these compressions
params.strainmeasure = 'pepicelli'; % which elastic constitutive model

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

% store the converged values of C and area0 for the elastic problem
params.area0 = pi*2*params.w*(itervars.r)/itervars.C;
params.C = itervars.C;

% NOTE: at this stage the initial guesses for r, z, psi and p0 are taken
% from the solution without elasticity. 

% initialize the surface strains and stresses
% NOTE: in the solution procedure, the dlams and dlamp (Newton-Raphson 
% update variables) are required to be equal at s=0. For the correct 
% solution, the initial guess must have equal lams and lamp at s=0
itervars.lamp = ones(params.N,1); itervars.lams = itervars.lamp;
itervars.sigmas = params.sigma*ones(params.N,1); 
itervars.sigmap = itervars.sigmas;

% store the coordinates of the reference shape
itervars.r_star = itervars.r; itervars.z_star = itervars.z;

% solve the elastic Young-Laplace equation for the given parameters
for ii = 1:length(params.fracm)

    % get the current value of the compression
    params.frac = params.fracm(ii);

    % solve the elastic Young-Laplace equation
    [itervars,params] = solve_young_laplace_elastic(itervars,params);

    % calculate the volume and the area
    wdef = params.w.*itervars.lams'/params.C; 
    volume = pi*wdef*(itervars.r.^2.*sin(itervars.psi));
    area = pi*2*wdef*(itervars.r);

    disp(['volume = ', num2str(volume,15)]);
    disp(['area = ', num2str(area,15)]);
    disp(['pressure = ', num2str(itervars.p0,15)]);

    % interpolate the numerical solutions on a finer grid. 
    % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
    % though all the points (see book of Trefethen on Spectral Methods in 
    % Matlab, page  63). For plotting purposes we use a simpler interpolation 
    % NOTE2: the interpolation is performed on the "numerical grid", 
    % this is not the actual value of s
    ss = linspace(params.s(1),params.s(end),params.Nplot)';
    rr = interp1(params.s,itervars.r,ss,'pchip');
    zz = interp1(params.s,itervars.z,ss,'pchip');

    % plot the droplet shape
    plot(rr,zz); 
    rmax = max([itervars.r_star',rr']);
    zmin = min([itervars.z_star',zz']);

    % rescale the plot
    xlim([0 1.2*rmax]);
    ylim([1.2*zmin 0]);
    set(gca,'DataAspectRatio',[1 1 1])

    % plot the surface stresses
    figure;
    plot(params.sdef,itervars.sigmas,'LineWidth',2); hold on
    plot(params.sdef,itervars.sigmap,'LineWidth',2);
    xlabel('s','FontSize',32);
    ylabel('\sigma','FontSize',32);
    legend('\sigma_s','\sigma_\phi','FontSize',24,'Location','northwest');
    xlim([0,params.sdef(end)])
    ax = gca; 
    ax.FontSize = 24;

end