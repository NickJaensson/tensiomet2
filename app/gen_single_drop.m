
close all; clear

% physical parameters
params_phys.sigma = 4;        % surface tension
params_phys.grav = 1.2;       % gravitational acceleration
params_phys.rneedle = 1.4;    % radius of the needle
params_phys.volume0 = 16;     % prescribed volume
params_phys.deltarho = 1.1;   % density difference

% numerical parameters
params_num.N = 40;          % resolution of the discretization for calculation
params_num.Nplot = 80;      % resolution of the discretization for plotting
params_num.eps_fw = 1e-12;  % convergence critertion forward: rms(u) < eps
params_num.maxiter = 100;   % maximum number of iteration steps

% calculate the Worthinton number
params_phys.Wo = params_phys.deltarho*params_phys.grav*params_phys.volume0/...
                                    (2*pi*params_phys.sigma*params_phys.rneedle);

% find an initial guess of the shape
shape_guess = guess_shape(params_phys,1000);

% get the differentation/integration matrices and the grid
vars_num = numerical_grid(params_num,[0,shape_guess.s(end)]);

% solve the Young-Laplace equation for the given parameters
vars_sol = solve_forward_young_laplace(params_phys, params_num, ...
    shape_guess, vars_num);

% update the values in the numerical grid
vars_num = update_numerical_grid(vars_sol,vars_num);

% calculate the volume and the area
volume = pi*vars_num.ws*(vars_sol.r.^2.*sin(vars_sol.psi));
area = pi*2*vars_num.ws*(vars_sol.r);

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

% determine the curvatures
% NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
% solved here by taking kappap(0) = kappas(0)
kappas = vars_num.Ds*vars_sol.psi;
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