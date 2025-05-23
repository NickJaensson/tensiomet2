% physical parameters for the simple droplet problem
params_phys.sigma = 20e-3;      % surface tension
params_phys.grav = 9.91;     % gravitational acceleration
params_phys.rneedle = 50e-3;  % radius of the needle
params_phys.volume0 = 0;   % prescribed volume
params_phys.deltarho = -1e3; % density difference
% params_phys.a = 43.47;
% params_phys.b = 0.006;
params_phys.a = -4.49;
params_phys.b = -729.56;
params_phys.pmax = 43.47;
params_phys.pmin = -0.002;

params_phys.impose_contact_angle = 1; % 0: set z=0 at wall
                                      % 1: set psi=contact_angle at wall
params_phys.contact_angle = 0;  

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);