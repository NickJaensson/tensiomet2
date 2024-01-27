% physical parameters for the simple droplet problem
params_phys.sigma = 4;      % surface tension
params_phys.grav = 1.2;     % gravitational acceleration
params_phys.rneedle = 1.4;  % radius of the needle
params_phys.volume0 = 16;   % prescribed volume
params_phys.deltarho = 1.1; % density difference

% numerical parameters
params_num.N = 40;          % grid points for calculation
params_num.eps_fw = 1e-12;  % convergence criterion forward: rms(u) < eps
params_num.maxiter = 100;   % maximum number of iteration steps

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);

% physical parameters for the elastic problem
params_phys.Kmod = 3;          % elastic dilational modulus
params_phys.Gmod = 2;          % elastic shear modulus
params_phys.compresstype = 1;  % 1: compress the volume other: compress the area
params_phys.fracm = [0.8];     % compute elastic stresses for these compressions
params_phys.strainmeasure = 'pepicelli'; % which elastic constitutive model