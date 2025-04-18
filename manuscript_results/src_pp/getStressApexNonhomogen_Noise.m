function [ApexStressTable] = getStressApexNonhomogen_Noise(DataStruc, NumSimuls, VarList)

ApexStress_fwd  = zeros(numel(VarList), 1);
ApexStress_TMet = zeros(numel(VarList), 1);
ApexStress_Circ = zeros(numel(VarList), 1);
StressParam_TMet = zeros(numel(VarList), 1);
StressParam_Circ = zeros(numel(VarList), 1);
K_vals = zeros(numel(VarList), 1);
G_vals = zeros(numel(VarList), 1);

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        K_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        G_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Gmod_dimal;
        % Get dilatational stress at the apex (sigmas or sigmap)
        ApexStress_fwd(idx,1)  = DataStruc.(VarList{i,j}).DefSt_fwd.sigmas_dimal(1,1);
        ApexStress_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Avg_dimal;
        ApexStressSTD_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Std_dimal;
        ApexStress_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.SigmaCircleFit_dimal_Avg;
        ApexStressSTD_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.SigmaCircleFit_dimal_Std;
        ApexStress_fwd_mat(i,j) = ApexStress_fwd(idx,1);
        ApexStress_TMet_mat(i,j) = ApexStress_TMet(idx,1);
        ApexStress_TMetSTD_mat(i,j) = ApexStressSTD_TMet(idx,1);
        ApexStress_Circ_mat(i,j) = ApexStress_Circ(idx,1);
        ApexStress_CircSTD_mat(i,j) = ApexStressSTD_Circ(idx,1);
        % strain param = (StressModel - StressReal) / StressReal
        StressParam_TMet(idx,1) = (ApexStress_TMet(idx,1) - ApexStress_fwd(idx,1))/ApexStress_fwd(idx,1);
        StressParam_Circ(idx,1) = (ApexStress_Circ(idx,1) - ApexStress_fwd(idx,1))/ApexStress_fwd(idx,1);
        StressParam_TMet_mat(i,j) = StressParam_TMet(idx,1);
        StressParam_Circ_mat(i,j) = StressParam_Circ(idx,1);
        % STD param = relative error = STD/avg
        StressParamSTD_TMet(idx,1) = ApexStressSTD_TMet(idx,1)/ApexStress_TMet(idx,1);
        StressParamSTD_Circ(idx,1) = ApexStressSTD_Circ(idx,1)/ApexStress_Circ(idx,1);
        StressParamSTD_TMet_mat(i,j) = StressParamSTD_TMet(idx,1);
        StressParamSTD_Circ_mat(i,j) = StressParamSTD_Circ(idx,1);
        idx = idx + 1;
    end
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
frac = DataStruc.(VarList{1,1}).Params_phys.frac;
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        ApexStress_fwd_mat2(i,j) = (ApexStress_fwd_mat(i,j) + ApexStress_fwd_mat(i+1,j) + ApexStress_fwd_mat(i,j+1) + ApexStress_fwd_mat(i+1,j+1))/4;
        ApexStress_TMet_mat2(i,j) = (ApexStress_TMet_mat(i,j) + ApexStress_TMet_mat(i+1,j) + ApexStress_TMet_mat(i,j+1) + ApexStress_TMet_mat(i+1,j+1))/4;
        ApexStress_Circ_mat2(i,j) = (ApexStress_Circ_mat(i,j) + ApexStress_Circ_mat(i+1,j) + ApexStress_Circ_mat(i,j+1) + ApexStress_Circ_mat(i+1,j+1))/4;
        StressParam_TMet_mat2(i,j) = (StressParam_TMet_mat(i,j) + StressParam_TMet_mat(i+1,j) + StressParam_TMet_mat(i,j+1) + StressParam_TMet_mat(i+1,j+1))/4;
        StressParam_Circ_mat2(i,j) = (StressParam_Circ_mat(i,j) + StressParam_Circ_mat(i+1,j) + StressParam_Circ_mat(i,j+1) + StressParam_Circ_mat(i+1,j+1))/4;
    end
end
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        ApexStress_TMetSTD_mat2(i,j) = (ApexStress_TMetSTD_mat(i,j) + ApexStress_TMetSTD_mat(i+1,j) + ApexStress_TMetSTD_mat(i,j+1) + ApexStress_TMetSTD_mat(i+1,j+1))/4;
        ApexStress_CircSTD_mat2(i,j) = (ApexStress_CircSTD_mat(i,j) + ApexStress_CircSTD_mat(i+1,j) + ApexStress_CircSTD_mat(i,j+1) + ApexStress_CircSTD_mat(i+1,j+1))/4;
        StressParamSTD_TMet_mat2(i,j) = (StressParamSTD_TMet_mat(i,j) + StressParamSTD_TMet_mat(i+1,j) + StressParamSTD_TMet_mat(i,j+1) + StressParamSTD_TMet_mat(i+1,j+1))/4;
        StressParamSTD_Circ_mat2(i,j) = (StressParamSTD_Circ_mat(i,j) + StressParamSTD_Circ_mat(i+1,j) + StressParamSTD_Circ_mat(i,j+1) + StressParamSTD_Circ_mat(i+1,j+1))/4;
    end
end

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            ApexStress_fwd2(idx, 1) = ApexStress_fwd_mat2(i,j);
            ApexStress_TMet2(idx, 1) = ApexStress_TMet_mat2(i,j);
            ApexStress_Circ2(idx, 1) = ApexStress_Circ_mat2(i,j);
            StressParam_TMet2(idx, 1) = StressParam_TMet_mat2(i,j);
            StressParam_Circ2(idx, 1) = StressParam_Circ_mat2(i,j);
        else
            ApexStress_fwd2(idx, 1) = ApexStress_fwd_mat(i,j);
            ApexStress_TMet2(idx, 1) = ApexStress_TMet_mat(i,j);
            ApexStress_Circ2(idx, 1) = ApexStress_Circ_mat(i,j);
            StressParam_TMet2(idx, 1) = StressParam_TMet_mat(i,j);
            StressParam_Circ2(idx, 1) = StressParam_Circ_mat(i,j);
        end
        idx = idx + 1;
    end
end
idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            ApexStressSTD_TMet2(idx, 1) = ApexStress_TMetSTD_mat2(i,j);
            ApexStress_Circ2(idx, 1) = ApexStress_CircSTD_mat2(i,j);
            StressParamSTD_TMet2(idx, 1) = StressParamSTD_TMet_mat2(i,j);
            StressParamSTD_Circ2(idx, 1) = StressParamSTD_Circ_mat2(i,j);
        else
            ApexStress_TMet2(idx, 1) = ApexStress_TMetSTD_mat(i,j);
            ApexStress_Circ2(idx, 1) = ApexStress_CircSTD_mat(i,j);
            StressParamSTD_TMet2(idx, 1) = StressParamSTD_TMet_mat(i,j);
            StressParamSTD_Circ2(idx, 1) = StressParamSTD_Circ_mat(i,j);
        end
        idx = idx + 1;
    end
end

%ApexStressTable = [G_vals, K_vals, ApexStress_fwd, ApexStress_TMet, ApexStress_Circ, StressParam_TMet, StressParam_Circ, StressParamSTD_TMet, StressParamSTD_Circ];
ApexStressTable = [G_vals, K_vals, ApexStress_fwd2, ApexStress_TMet2, ApexStress_Circ2, StressParam_TMet2, StressParam_Circ2, StressParamSTD_TMet2, StressParamSTD_Circ2];

end