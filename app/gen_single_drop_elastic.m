% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% physical parameters
params_phys.sigma = 4;        % surface tension
params_phys.grav = 1.2;       % gravitational acceleration
params_phys.rneedle = 1.4;    % radius of the needle
params_phys.volume0 = 16;     % prescribed volume
params_phys.deltarho = 1.1;   % density difference


% physical parameters for the elastic problem
params_phys.Kmod = 3;          % elastic dilational modulus
params_phys.Gmod = 2;          % elastic shear modulus
params_phys.compresstype = 1;  % 1: compress the volume other: compress the area
params_phys.fracm = [0.8];     % compute elastic stresses for these compressions
params_phys.strainmeasure = 'pepicelli'; % which elastic constitutive model

% numerical parameters
params_num.N = 40;          % resolution of the discretization for calculation
params_num.Nplot = 80;      % resolution of the discretization for plotting
params_num.maxiter = 100;   % maximum number of iteration steps
params_num.eps_fw = 1e-12;  % convergence critertion forward: rms(u) < eps

% calculate the Worthinton number
params_phys.Wo = params_phys.deltarho*params_phys.grav*params_phys.volume0/...
                                    (2*pi*params_phys.sigma*params_phys.rneedle);

% solve the Young-Laplace equation for the given parameters
[vars_sol,vars_num] = solve_forward_young_laplace(params_phys, params_num);

% calculate the volume and the area
volume = pi*vars_num.w*(vars_sol.r.^2.*sin(vars_sol.psi))/vars_sol.C;
area = pi*2*vars_num.w*(vars_sol.r)/vars_sol.C;

disp(['volume = ', num2str(volume,15)]);
disp(['area = ', num2str(area,15)]);
disp(['pressure = ', num2str(vars_sol.p0,15)]);

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(vars_num.s(1),vars_num.s(end),params_num.Nplot)';
rr = interp1(vars_num.s,vars_sol.r,ss,'pchip');
zz = interp1(vars_num.s,vars_sol.z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rr',zz','b');
plot(rr',zz','b');
set(gca,'DataAspectRatio',[1 1 1])

% store the converged values of C and area0 for the elastic problem
params_phys.area0 = pi*2*vars_num.w*(vars_sol.r)/vars_sol.C;
vars_num.C = vars_sol.C;

% NOTE: at this stage the initial guesses for r, z, psi and p0 are taken
% from the solution without elasticity. 

% initialize the surface strains and stresses
% NOTE: in the solution procedure, the dlams and dlamp (Newton-Raphson 
% update variables) are required to be equal at s=0. For the correct 
% solution, the initial guess must have equal lams and lamp at s=0
vars_sol.lamp = ones(params_num.N,1); vars_sol.lams = vars_sol.lamp;
vars_sol.sigmas = params_phys.sigma*ones(params_num.N,1); 
vars_sol.sigmap = vars_sol.sigmas;

% store the coordinates of the reference shape
vars_sol.r_star = vars_sol.r; vars_sol.z_star = vars_sol.z;

% solve the elastic Young-Laplace equation for the given parameters
for ii = 1:length(params_phys.fracm)

    % get the current value of the compression
    params_phys.frac = params_phys.fracm(ii);

    % solve the elastic Young-Laplace equation
    [vars_sol,vars_num] = solve_forward_young_laplace_elastic(vars_sol, params_phys, params_num, vars_num);

    % calculate the volume and the area
    volume = pi*vars_num.wdef*(vars_sol.r.^2.*sin(vars_sol.psi));
    area = pi*2*vars_num.wdef*(vars_sol.r);

    disp(['volume = ', num2str(volume,15)]);
    disp(['area = ', num2str(area,15)]);
    disp(['pressure = ', num2str(vars_sol.p0,15)]);

    % interpolate the numerical solutions on a finer grid. 
    % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
    % though all the points (see book of Trefethen on Spectral Methods in 
    % Matlab, page  63). For plotting purposes we use a simpler interpolation 
    % NOTE2: the interpolation is performed on the "numerical grid", 
    % this is not the actual value of s
    ss = linspace(vars_num.s(1),vars_num.s(end),params_num.Nplot)';
    rr = interp1(vars_num.s,vars_sol.r,ss,'pchip');
    zz = interp1(vars_num.s,vars_sol.z,ss,'pchip');

    % plot the droplet shape
    plot(rr,zz); 
    rmax = max([vars_sol.r_star',rr']);
    zmin = min([vars_sol.z_star',zz']);

    % rescale the plot
    xlim([0 1.2*rmax]);
    ylim([1.2*zmin 0]);
    set(gca,'DataAspectRatio',[1 1 1])

    % plot the surface stresses
    figure;
    plot(vars_num.sdef,vars_sol.sigmas,'LineWidth',2); hold on
    plot(vars_num.sdef,vars_sol.sigmap,'LineWidth',2);
    xlabel('s','FontSize',32);
    ylabel('\sigma','FontSize',32);
    legend('\sigma_s','\sigma_\phi','FontSize',24,'Location','northwest');
    xlim([0,vars_num.sdef(end)])
    ax = gca; ax.FontSize = 24;

    % determine the curvatures
    % NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
    % solved here by taking kappap(0) = kappas(0)
    kappas = vars_num.Ddef*vars_sol.psi;
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

end