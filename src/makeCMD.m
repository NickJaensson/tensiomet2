function [ taus, taup ] = makeCMD(params_phys, psi, r, z, vars_num, p0)
% compute stresses from the shape and a pressure following Danov et al. 
% Advances in Colloid and Interface Science 233 (2016) 223?239
        
    fac = params_phys.deltarho*params_phys.grav;

    Vi = pi*vars_num.wsmat*(r.^2.*sin(psi));

    L = pi*2*r.*sin(psi);
    R = pi*(r.^2*p0 - fac*z.*r.^2);

    taus = zeros(vars_num.N,1);
    taus(2:end) = ( R(2:end) + fac*Vi(2:end) ) ./ L(2:end);

    % regularize at the apex with dsigma^s(0)/ds =0
    taus(1) = -(vars_num.Ds(1,2:end)*taus(2:end))./vars_num.Ds(1,1);

    taup = (p0-fac*z-taus.*(vars_num.Ds*psi)).*r./sin(psi);
    taup(1) = taus(1);

end

