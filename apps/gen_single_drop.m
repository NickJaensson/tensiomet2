
close all; clear

addpath('../src/')

% physical parameters
params.sigma = 72;        % surface tension [mN/m]
params.grav = 9.807e3;    % gravitational acceleration [mm/s^2]
params.rneedle = 1;       % radius of the needle [mm]
params.volume0 = 32;       % prescribed volume in mm^3
params.deltarho = 1e-3;   % density difference [10^6 kg/m^3]

% numerical parameters
params.N = 40;          % resolution of the discretization for calculation
params.Nplot = 80;      % resolution of the discretization for plotting
params.Ncheb = 10;      % number of Chebyshev to describe the shape
params.alpha = 1;       % relaxation parameter in the Newton-Raphson scheme

tic

% NOTE: the calculation is done in dimensionless form, using the 
% dimensionless surface tension sigma' and volume V'

% calculate the dimensionless quantities
params.sigmaprime = params.sigma/ ...
                            (params.deltarho*params.grav*params.rneedle^2);
params.volume0prime = params.volume0/params.rneedle^3;

% find an initial guess of the shape
[r_guess, z_guess, s_guess] = guess_shape(params,1000);

smax = s_guess(end); % the total length of the 1D domain

% get the differentation/integration matrices and the grid
[params.D,~,params.w,s] = dif1D('cheb',0,smax,params.N,5);

% interpolate the shape in the Chebyshev points
r = interp1(s_guess,r_guess,s);
z = interp1(s_guess,z_guess,s);

psi = atan2(params.D*z,params.D*r);      % intial psi value 
C = 1;                     % initial stretch parameter
p0 = 2*params.rneedle*params.sigmaprime;   % predict the pressure

u = ones(3*params.N+2,1);

iter = 0; crash = 0; 

while rms(u) > 1e-10

    iter = iter + 1;
    
    if iter > 1200 
        warning('iter > 12000!');
        crash = 1; break;
    end
    
    itervars.r = r; itervars.z = z; itervars.psi = psi;
    itervars.C = C; itervars.p0 = p0;  

    [A,b] = jacobian_rhs_simple(params,itervars);
    
    % solve the system of equations
    u = A\b;
    
    % update variables
    r   = r   + params.alpha*u(1:params.N);
    z   = z   + params.alpha*u(params.N+1:2*params.N);
    psi = psi + params.alpha*u(2*params.N+1:3*params.N); 
    C   = C   + params.alpha*u(3*params.N+1);
    p0  = p0  + params.alpha*u(3*params.N+2);
    
    if rms(b) > 1e3
        crash = 1; break;
    end
    
    fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

end

toc

% compute volume and area (scaled back to dimensionfull)
disp(['volume = ', num2str(params.rneedle^3*pi*params.w*(r.^2.*sin(psi))/C,15),' mm^3']);
disp(['area = ', num2str(params.rneedle^2*pi*2*params.w*(r)/C,15),' mm^2']);
disp(['pressure = ', num2str(params.deltarho*params.grav*params.rneedle*p0,15),' Pa']);

% % plot the shape of the drop on the numerical grid
% figure; hold on
% scatter(rneedle*r',rneedle*z','b');
% plot(rneedle*r',rneedle*z','b');
% set(gca,'DataAspectRatio',[1 1 1])

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(s(1),s(end),params.Nplot)';
rr = interp1(s,r,ss,'pchip');
zz = interp1(s,z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(params.rneedle*rr',params.rneedle*zz','b');
plot(params.rneedle*rr',params.rneedle*zz','b');
set(gca,'DataAspectRatio',[1 1 1])