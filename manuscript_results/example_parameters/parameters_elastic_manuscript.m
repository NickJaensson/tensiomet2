% physical parameters for the elastic droplet problem

% dimensionfull input parameters
Kmod = 101;        % elastic dilational modulus [mN/m]
Gmod = 2;        % elastic shear modulus [mN/m]
compresstype = 2;  % 1: compress the volume other: compress the area
frac  = 0.8;       % compute elastic stresses for this compression
strainmeasure = 'pepicelli'; % which elastic constitutive model (linear_hookean, hencky, balemans, or pepicelli) 

% dimensional input parameters, for saving
params_phys.Kmod_dimal = Kmod;
params_phys.Gmod_dimal = Gmod;

% dimensionless input parameters for calculation
params_phys.Kmod = Kmod/(deltarho*abs(grav)*rneedle^2);
params_phys.Gmod = Gmod/(deltarho*abs(grav)*rneedle^2);
params_phys.compresstype = compresstype;  
params_phys.frac  = frac;
params_phys.strainmeasure = strainmeasure; 

