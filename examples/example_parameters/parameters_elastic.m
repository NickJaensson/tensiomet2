% physical parameters for the elastic problem
params_phys.Kmod = 3;          % elastic dilational modulus
params_phys.Gmod = 2;          % elastic shear modulus
params_phys.compresstype = 1;  % 1: compress the volume other: compress the area
params_phys.frac  = 0.8;       % compute elastic stresses for this compression
params_phys.strainmeasure = 'pepicelli'; % which elastic constitutive model