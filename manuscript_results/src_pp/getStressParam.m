function [StressParam_table] = getStressParam(DataStruc, VarList)

stress_param = zeros(numel(VarList), 1);
K_vals = zeros(numel(VarList), 1);
G_vals = zeros(numel(VarList), 1);

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        stress_s = DataStruc.(VarList{i,j}).DefSt_fwd.sigmas_dimal;
        stress_p = DataStruc.(VarList{i,j}).DefSt_fwd.sigmap_dimal;
        stress_param(idx, 1) = (stress_p(end) - stress_s(end))/(stress_s(1));
        stress_param_mat(i,j) = stress_param(idx, 1);
        K_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        G_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Gmod_dimal;
        idx = idx + 1;
    end
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
frac = DataStruc.(VarList{1,1}).Params_phys.frac;
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        if  Sigma == 20 && frac == 0.8  %excption for inaccessible geometries, average less cells
            if i == size(VarList,1)-1 
                stress_param_mat2(i,j) = stress_param_mat(i,j);
            else
                stress_param_mat2(i,j) = (stress_param_mat(i,j) + stress_param_mat(i+1,j) + stress_param_mat(i,j+1) + stress_param_mat(i+1,j+1))/4;
            end
        else
            stress_param_mat2(i,j) = (stress_param_mat(i,j) + stress_param_mat(i+1,j) + stress_param_mat(i,j+1) + stress_param_mat(i+1,j+1))/4;
        end
    end
end

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            stress_param2(idx, 1) = stress_param_mat2(i,j);
        else
            stress_param2(idx, 1) = stress_param_mat(i,j);
        end
        idx = idx + 1;
    end
end

%StressParam_table = [G_vals, K_vals, stress_param];
StressParam_table = [G_vals, K_vals, stress_param2];


end