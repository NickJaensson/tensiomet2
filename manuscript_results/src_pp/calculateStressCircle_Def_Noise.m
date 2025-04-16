function [DataStruc] = calculateStressCircle_Def_Noise(DataStruc, VarList, NumSimuls, fraction, show)

R_circlefit = zeros(1,NumSimuls);
sigma_circlefit = zeros(1,NumSimuls);
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        for kk = 1:NumSimuls

            % Case 1: subhemispherical bubble
            if min(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal) == DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal(1,1)
                Lfrac = ceil((1-fraction)*length(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal)/2);
                [r_circlefit, z_circlefit, R_circlefit(kk)] = getCircle(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal, DataStruc.(VarList{i,j}).DefSt_inv{kk}.z_noise_dimal, Lfrac);
                % Case 2: suphemispherical bubble
            else
                idx_min = find(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal == min(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal));
                idx_max = find(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal == max(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal));
                Lfrac = ceil((1-fraction)*length(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal(idx_min:idx_max))/2);
                if rem(length(DataStruc.(VarList{i,j}).RefSt_inv{kk}.r_noise_dimal(idx_min:idx_max)), 2) == 0
                    idx_max = idx_max - 1;
                end
                [r_circlefit, z_circlefit, R_circlefit(kk)] = getCircle(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal(idx_min:idx_max), ...
                    DataStruc.(VarList{i,j}).DefSt_inv{kk}.z_noise_dimal(idx_min:idx_max), Lfrac);
            end

            % Calculate surface stress using circle fitting.
            % Use prescribed density difference for air-water
            grav = 9.807;    % gravitational acceleration [m/s^2]
            deltarho = 1000;   % density difference [kg/m^3]
            sigma_circlefit(kk) = (DataStruc.(VarList{i,j}).DefSt_inv{kk}.p0_dimal - deltarho*grav*min(DataStruc.(VarList{i,j}).DefSt_inv{kk}.z_noise_dimal)/1000)*R_circlefit(kk)/2;

            % save data
            DataStruc.(VarList{i,j}).DefSt_inv{kk}.RadiusCircleFit_dimal = R_circlefit(kk);
            DataStruc.(VarList{i,j}).DefSt_inv{kk}.SigmaCircleFit_dimal = sigma_circlefit(kk);

            % plot (can comment out)
            if show == 1
                figure(100)
                plot(DataStruc.(VarList{i,j}).DefSt_inv{kk}.r_noise_dimal, DataStruc.(VarList{i,j}).DefSt_inv{kk}.z_noise_dimal, 'bo'); hold on
                plot(r_circlefit, z_circlefit, 'r-');
            end
        end

        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.RadiusCircleFit_dimal_Avg = mean(R_circlefit);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.RadiusCircleFit_dimal_Std = std(R_circlefit);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.SigmaCircleFit_dimal_Avg = mean(sigma_circlefit);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.SigmaCircleFit_dimal_Std = std(sigma_circlefit);
    end
end

end