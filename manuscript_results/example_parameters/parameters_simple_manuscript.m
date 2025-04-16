% physical parameters for the simple droplet problem

% dimensionfull input parameters
Wo = 0.1;
Ar = 5;
sigma = 60;        % surface tension [mN/m]
grav = 9.807e3;    % gravitational acceleration [mm/s^2]
deltarho = 1e-3;   % density difference [10^6 kg/m^3]

% calculation of dimensionfull parameters
rneedle = sqrt((Wo*sigma/(deltarho*abs(grav)*Ar))); % radius of the needle [mm]
volume0 = Ar*rneedle^3;  % prescribed volume in mm^3

% dimensional input parameters, for saving
params_phys.sigma_dimal = sigma;
params_phys.grav_dimal = grav;
params_phys.rneedle_dimal = rneedle;
params_phys.volume0_dimal = volume0;
params_phys.deltarho_dimal = deltarho;

% dimensionless input parameters for calculation
params_phys.sigma = sigma/(deltarho*abs(grav)*rneedle^2);
params_phys.grav = grav/abs(grav);
params_phys.rneedle = rneedle/rneedle;
params_phys.volume0 = volume0/rneedle^3;
params_phys.deltarho = deltarho/deltarho;

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);

% Dimensionless numbers for the paper
params_phys.Wo_paper = Wo;
params_phys.Ar_paper = Ar;



