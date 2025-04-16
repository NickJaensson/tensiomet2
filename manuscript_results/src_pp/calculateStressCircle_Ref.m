function [DataStruc] = calculateStressCircle_Ref(DataStruc, VarList, fraction, show)

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)

        % Case 1: subhemispherical bubble
        if min(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal) == DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal(1,1)
            Lfrac = ceil((1-fraction)*length(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal)/2);
            [r_circlefit, z_circlefit, R_circlefit] = getCircle(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal, DataStruc.(VarList{i,j}).RefSt_inv.z_noise_dimal, Lfrac);
        % Case 2: suphemispherical bubble
        else
            idx_min = find(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal == min(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal));
            idx_max = find(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal == max(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal));
            Lfrac = ceil((1-fraction)*length(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal(idx_min:idx_max))/2);
            [r_circlefit, z_circlefit, R_circlefit] = getCircle(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal(idx_min:idx_max), ...
                DataStruc.(VarList{i,j}).RefSt_inv.z_noise_dimal(idx_min:idx_max), Lfrac);
        end

        % Calculate surface stress using circle fitting.
        % Use prescribed density difference for air-water
        grav = 9.807;    % gravitational acceleration [m/s^2]
        deltarho = 1000;   % density difference [kg/m^3]
        % p0 in Pa, deltarho*grav in SI units, z in mm (divide by 1000 to convert to m). R_circlefit in mm, so sigma_circlefit in mN/m
        sigma_circlefit = (DataStruc.(VarList{i,j}).RefSt_inv.p0_dimal - deltarho*grav*min(DataStruc.(VarList{i,j}).RefSt_inv.z_noise_dimal)/1000)*R_circlefit/2;

        % save data
        DataStruc.(VarList{i,j}).RefSt_inv.RadiusCircleFit_dimal = R_circlefit;
        DataStruc.(VarList{i,j}).RefSt_inv.SigmaCircleFit_dimal = sigma_circlefit;

        % plot (can comment out)
        if show == 1
            figure(100)
            plot(DataStruc.(VarList{i,j}).RefSt_inv.r_noise_dimal, DataStruc.(VarList{i,j}).RefSt_inv.z_noise_dimal, 'bo'); hold on
            plot(r_circlefit, z_circlefit, 'r-');
        end

    end
end

end