
close all; clear

addpath('../src/')

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

% calculate the Worthinton number
params.Wo = params.deltarho*params.grav*params.volume0/...
                                        (2*pi*params.sigma*params.rneedle);

disp(['Wo = ', num2str(params.Wo,3)]);

% find an initial guess of the shape
[r_guess, z_guess, s_guess] = guess_shape(params,1000);

smax = s_guess(end); % the total length of the 1D domain

% get the differentation/integration matrices and the grid
[params.D,~,params.w,s] = dif1D('cheb',0,smax,params.N,5);

% interpolate the shape in the Chebyshev points
r = interp1(s_guess,r_guess,s);
z = interp1(s_guess,z_guess,s);

psi = atan2(params.D*z,params.D*r);   % intial psi value 
C = 1;                                % initial stretch parameter
p0 = 2*params.sigma/params.rneedle;   % initial pressure
u = ones(3*params.N+2,1);             % initial solution vector

% store some variables for the iteration
iter = 0;
itervars.r = r; itervars.z = z; itervars.psi = psi;
itervars.C = C; itervars.p0 = p0; 

% start the Newton-Raphson iteration
while rms(u) > params.eps

    iter = iter + 1;
    
    if iter > params.maxiter
        error('Iteration did not converge!')
    end    

    % build the Jacobian and RHS
    [A,b] = jacobian_rhs_simple(params,itervars);
    
    % solve the system of equations
    u = A\b;
    
    % update variables
    itervars.r   = itervars.r   + params.alpha*u(1:params.N);
    itervars.z   = itervars.z   + params.alpha*u(params.N+1:2*params.N);
    itervars.psi = itervars.psi + params.alpha*u(2*params.N+1:3*params.N); 
    itervars.C   = itervars.C   + params.alpha*u(3*params.N+1);
    itervars.p0  = itervars.p0  + params.alpha*u(3*params.N+2);    
    
    fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

end

% extract the solution variables
r = itervars.r; z = itervars.z; psi = itervars.psi;
C = itervars.C; p0 = itervars.p0;

% calculate the volume and the area
volume = pi*params.w*(r.^2.*sin(psi))/C;
area = pi*2*params.w*(r)/C;

disp(['volume = ', num2str(volume,15)]);
disp(['area = ', num2str(area,15)]);
disp(['pressure = ', num2str(p0,15)]);

% store the converged values of C and area0 for the elastic problem
params.area0 = pi*2*params.w*(r)/C;
params.C = C;

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(s(1),s(end),params.Nplot)';
rr = interp1(s,r,ss,'pchip');
zz = interp1(s,z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rr',zz','b');
plot(rr',zz','b');
set(gca,'DataAspectRatio',[1 1 1])