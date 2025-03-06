function [volume, area] = ...
    calculate_volume_area(vars_sol, vars_num, verbose)
    % CALCULATE_VOLUME_AREA Computes volume and surface area.
    %
    % INPUTS:
    %   vars_sol  - Structure with solution variables
    %   vars_num  - Structure with numerical variables
    %   verbose   - Boolean flag to display results if true
    %
    % OUTPUTS:
    %   volume    - Computed volume
    %   area      - Computed surface area

    % calculate the volume and the area
    volume = pi*vars_num.ws*(vars_sol.r.^2.*sin(vars_sol.psi));
    area = pi*2*vars_num.ws*(vars_sol.r);

    if verbose
        disp(['volume = ', num2str(volume,15)]);
        disp(['area = ', num2str(area,15)]);
        disp(['pressure = ', num2str(vars_sol.p0,15)]);
    end

end
