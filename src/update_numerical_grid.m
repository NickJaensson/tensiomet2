function [vars_num] = update_numerical_grid(vars_sol,vars_num)

    % the integration and differentation matrices in the solution state
    vars_num.ws = vars_num.w/vars_sol.C; 
    vars_num.Ds = vars_sol.C*vars_num.D; 
    vars_num.s = vars_num.s0/vars_sol.C;

end