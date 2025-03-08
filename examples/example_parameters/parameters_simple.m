% physical parameters for the simple droplet problem
params_phys.sigma = 4;      % surface tension
params_phys.grav = 1.2;     % gravitational acceleration
params_phys.rneedle = 1.4;  % radius of the needle
params_phys.volume0 = 16;   % prescribed volume
params_phys.deltarho = 1.1; % density difference

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);