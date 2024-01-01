function [ss,rr,zz] = interpolate_solutions(vars_sol, vars_num, params_num)

    % interpolate the numerical solutions on a uniform grid.
    % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
    % though all the points (see book of Trefethen on Spectral Methods in 
    % Matlab, page  63). For plotting purposes we use a simpler interpolation 
    ss = linspace(vars_num.s(1),vars_num.s(end),params_num.Nplot)';
    rr = interp1(vars_num.s,vars_sol.r,ss,'pchip');
    zz = interp1(vars_num.s,vars_sol.z,ss,'pchip');

end

