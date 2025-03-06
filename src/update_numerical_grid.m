function [vars_num] = update_numerical_grid(vars_sol, vars_num, elastic)
    % UPDATE_NUMERICAL_GRID Updates numerical grid and matrices based on 
    % solution state.
    %
    % INPUTS:
    %   vars_sol  - Solution variables
    %   vars_num  - Numerical variables
    %   elastic   - Boolean flag indicating whether the system is elastic.
    %
    % OUTPUTS:
    %   vars_num  - Updated numerical variables

    % the integration and differentation matrices in the solution state
    if ~elastic
        vars_num.ws = vars_num.w/vars_sol.C; 
        vars_num.Ds = vars_sol.C*vars_num.D; 
        vars_num.s = vars_num.s0/vars_sol.C;
        vars_num.wsmat = vars_num.wmat/vars_sol.C;
        vars_num.C = vars_sol.C;
    else
        % the integration and differentation matrices in the deformed state
        % NOTE: this construction of Ds is simlar to first applying D*f/C,
        % and then dividing the components by the components of (1/lams)
        vars_num.ws = vars_num.w.*vars_sol.lams'/vars_num.C; 
        vars_num.Ds = vars_num.C*vars_num.D.*repelem((1./vars_sol.lams),1,vars_num.N); 
    
        % get the integration matrix in th reference state
        % NOTE: this is NOT the integration matrix for the deformed state!
        vars_num.wsstarmat = vars_num.wmat/vars_num.C;
    
        % compute the value of s in the deformed state
        vars_num.s = vars_num.wsstarmat*vars_sol.lams;
    end    

end
