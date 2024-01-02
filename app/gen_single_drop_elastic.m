% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

% first calculate the reference state
gen_single_drop;

% physical parameters for the elastic problem
params_phys.Kmod = 3;          % elastic dilational modulus
params_phys.Gmod = 2;          % elastic shear modulus
params_phys.compresstype = 1;  % 1: compress the volume other: compress the area
params_phys.fracm = [0.8];     % compute elastic stresses for these compressions
params_phys.strainmeasure = 'pepicelli'; % which elastic constitutive model

% store the converged values of C and area0 for the elastic problem
params_phys.area0 = area;
vars_num.C = vars_sol.C;

% NOTE: at this stage the initial guesses for r, z, psi and p0 are taken
% from the solution without elasticity. 

% initialize the surface strains and stresses
% NOTE: in the solution procedure, the dlams and dlamp (Newton-Raphson 
% update variables) are required to be equal at s=0. For the correct 
% solution, the initial guess must have equal lams and lamp at s=0
vars_sol.lamp = ones(params_num.N,1); vars_sol.lams = vars_sol.lamp;
vars_sol.sigmas = params_phys.sigma*ones(params_num.N,1); 
vars_sol.sigmap = vars_sol.sigmas;

vars_sol.r_star = vars_sol.r; vars_sol.z_star = vars_sol.z;

for ii = 1:length(params_phys.fracm)

    params_phys.frac = params_phys.fracm(ii);

    [vars_sol,vars_num] = ...
        solve_forward_young_laplace_elastic(vars_sol, params_phys, ...
                                            params_num, vars_num);

    vars_num = update_numerical_grid(vars_sol, vars_num, true);

    [volume,area] = calculate_volume_area(vars_sol, vars_num, true);

    plot_shape(vars_sol.r, vars_sol.z, 1);

    plot_surface_stress(vars_num.s, vars_sol.sigmas, vars_sol.sigmap, 2);

    [kappas,kappap] = find_curvature(vars_sol, vars_num);

    plot_curvature(vars_sol.z, kappas, kappap, 3);

end