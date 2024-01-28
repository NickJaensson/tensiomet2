Vo = 5;  % called Ar in the manuscript
Bo = 1;  % called Wo in the manuscript

% physical parameters for the simple droplet problem
params_phys.sigma = 60;       % surface tension
params_phys.grav = 9.807e3;   % gravitational acceleration
params_phys.deltarho = 1e-3;  % density difference

params_phys.rneedle = sqrt(Bo * params_phys.sigma / ...
    (params_phys.deltarho * params_phys.grav * Vo));
params_phys.volume0 = Vo*params_phys.rneedle^3;

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);