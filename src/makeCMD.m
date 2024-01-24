function [ taus, taup ] = makeCMD(params_phys, vars_sol, vars_num)
% compute stresses from the shape and a pressure following Danov et al. 
% Advances in Colloid and Interface Science 233 (2016) 223?239

    psi = vars_sol.psi; 
    r = vars_sol.r;
    z = vars_sol.z;
    p0 = vars_sol.p0;

    fac = params_phys.deltarho*params_phys.grav;
    d = vars_num.Ds;

    % compute body force
    intH = [0; d(2:end,2:end)\(r(2:end).*r(2:end).*sin(psi(2:end)))];
    % line force for unknown surface stress
    L = 2*r.*sin(psi);
    % known forces
    rhs = fac*intH - fac*z.*r.^2 + r.^2*p0;
    
    taus = zeros(length(r),1);
    % force balance to obtain the meridional surface stress
    taus(2:end) = rhs(2:end)./L(2:end);
    % regularize at the apex with d sigma^s(0)/ds =0
    taus(1) = -(d(1,2:end)*taus(2:end))./d(1,1);
    
    % formulation like in Danov et al.
    taup = (p0-fac*z-taus.*(d*psi)).*r./sin(psi);
    taup(1) = taus(1);

end

