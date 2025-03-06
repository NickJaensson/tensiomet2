function [kappas, kappap] = find_curvature(vars_sol, vars_num)
    % FIND_CURVATURE Computes the curvature values.
    %
    % INPUTS:
    %   vars_sol  - Structure with solution variables
    %   vars_num  - Structure with numerical variables
    %
    % OUTPUTS:
    %   kappas    - Curvature along the s direction
    %   kappap    - Curvature along the phi direction

    % determine the curvatures
    % NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
    % solved here by taking kappap(0) = kappas(0)
    kappas = vars_num.Ds*vars_sol.psi;
    kappap = kappas;
    kappap(2:end) = sin(vars_sol.psi(2:end))./vars_sol.r(2:end);
    
end
