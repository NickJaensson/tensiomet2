function [volume,area] = calculate_volume_area(vars_sol,vars_num,print)

    % calculate the volume and the area
    volume = pi*vars_num.ws*(vars_sol.r.^2.*sin(vars_sol.psi));
    area = pi*2*vars_num.ws*(vars_sol.r);

    if print
        disp(['volume = ', num2str(volume,15)]);
        disp(['area = ', num2str(area,15)]);
        disp(['pressure = ', num2str(vars_sol.p0,15)]);
    end

end