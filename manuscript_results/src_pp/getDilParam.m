function [DilParam_table] = getDilParam(DataStruc, VarList)

dil_param = zeros(numel(VarList), 1);
K_vals = zeros(numel(VarList), 1);
G_vals = zeros(numel(VarList), 1);

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        dil = DataStruc.(VarList{i,j}).DefSt_fwd.lams.*DataStruc.(VarList{i,j}).DefSt_fwd.lamp;
        dil_param(idx, 1) = (max(dil) - min(dil))/(1-DataStruc.(VarList{i,j}).Params_phys.frac);
        dil_param_mat(i,j) = dil_param(idx, 1);
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
                dil_param_mat2(i,j) = dil_param_mat(i,j);
            else
                dil_param_mat2(i,j) = (dil_param_mat(i,j) + dil_param_mat(i+1,j) + dil_param_mat(i,j+1) + dil_param_mat(i+1,j+1))/4;
            end
        else
            dil_param_mat2(i,j) = (dil_param_mat(i,j) + dil_param_mat(i+1,j) + dil_param_mat(i,j+1) + dil_param_mat(i+1,j+1))/4;
        end
    end
end

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            dil_param2(idx, 1) = dil_param_mat2(i,j);
        else
            dil_param2(idx, 1) = dil_param_mat(i,j);
        end
        idx = idx + 1;
    end
end

%DilParam_table = [G_vals, K_vals, dil_param];
DilParam_table = [G_vals, K_vals, dil_param2];

end