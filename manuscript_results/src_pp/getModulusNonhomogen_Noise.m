function [ModulusTable] = getModulusNonhomogen_Noise(DataStruc, NumSimuls, VarList)

Kmod_fwd  = zeros(numel(VarList), 1);
Kmod_TMet = zeros(numel(VarList), 1);
KmodSTD_TMet = zeros(numel(VarList), 1);
Kmod_Circ = zeros(numel(VarList), 1);
KmodSTD_Circ = zeros(numel(VarList), 1);
KmodParam_TMet = zeros(numel(VarList), 1);
KmodParamSTD_TMet = zeros(numel(VarList), 1);
KmodParam_Circ = zeros(numel(VarList), 1);
KmodParamSTD_Circ = zeros(numel(VarList), 1);
K_vals = zeros(numel(VarList), 1);
G_vals = zeros(numel(VarList), 1);

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        % prescribed moduli, for forward problem
        K_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        G_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Gmod_dimal;
        % Get moduli
        Kmod_fwd(idx,1)  = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        Kmod_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Kmod_Avg_dimal;
        KmodSTD_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Kmod_Std_dimal;
        Kmod_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.KmodCircleFit_dimal_Avg;
        KmodSTD_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.KmodCircleFit_dimal_Std;
        Kmod_fwd_mat(i,j) = Kmod_fwd(idx,1);
        Kmod_TMet_mat(i,j) = Kmod_TMet(idx,1);
        Kmod_TMetSTD_mat(i,j) = KmodSTD_TMet(idx,1);
        Kmod_Circ_mat(i,j) = Kmod_Circ(idx,1);
        Kmod_CircSTD_mat(i,j) = KmodSTD_Circ(idx,1);
        % modulus param = (ModulusModel - ModulusReal) / ModulusReal
        KmodParam_TMet(idx,1) = (Kmod_TMet(idx,1) - Kmod_fwd(idx,1))/Kmod_fwd(idx,1);
        KmodParam_Circ(idx,1) = (Kmod_Circ(idx,1) - Kmod_fwd(idx,1))/Kmod_fwd(idx,1);
        KmodParam_TMet_mat(i,j) = KmodParam_TMet(idx,1);
        KmodParam_Circ_mat(i,j) = KmodParam_Circ(idx,1);
        % STD param = relative error = STD/avg
        KmodParamSTD_TMet(idx,1) = KmodSTD_TMet(idx,1)/Kmod_TMet(idx,1);
        KmodParamSTD_Circ(idx,1) = KmodSTD_Circ(idx,1)/Kmod_Circ(idx,1);
        KmodParamSTD_TMet_mat(i,j) = KmodParamSTD_TMet(idx,1);
        KmodParamSTD_Circ_mat(i,j) = KmodParamSTD_Circ(idx,1);
        idx = idx + 1;
    end
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
frac = DataStruc.(VarList{1,1}).Params_phys.frac;
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        Kmod_fwd_mat2(i,j) = (Kmod_fwd_mat(i,j) + Kmod_fwd_mat(i+1,j) + Kmod_fwd_mat(i,j+1) + Kmod_fwd_mat(i+1,j+1))/4;
        Kmod_TMet_mat2(i,j) = (Kmod_TMet_mat(i,j) + Kmod_TMet_mat(i+1,j) + Kmod_TMet_mat(i,j+1) + Kmod_TMet_mat(i+1,j+1))/4;
        Kmod_Circ_mat2(i,j) = (Kmod_Circ_mat(i,j) + Kmod_Circ_mat(i+1,j) + Kmod_Circ_mat(i,j+1) + Kmod_Circ_mat(i+1,j+1))/4;
        KmodParam_TMet_mat2(i,j) = (KmodParam_TMet_mat(i,j) + KmodParam_TMet_mat(i+1,j) + KmodParam_TMet_mat(i,j+1) + KmodParam_TMet_mat(i+1,j+1))/4;
        KmodParam_Circ_mat2(i,j) = (KmodParam_Circ_mat(i,j) + KmodParam_Circ_mat(i+1,j) + KmodParam_Circ_mat(i,j+1) + KmodParam_Circ_mat(i+1,j+1))/4;
    end
end
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        Kmod_TMetSTD_mat2(i,j) = (Kmod_TMetSTD_mat(i,j) + Kmod_TMetSTD_mat(i+1,j) + Kmod_TMetSTD_mat(i,j+1) + Kmod_TMetSTD_mat(i+1,j+1))/4;
        Kmod_CircSTD_mat2(i,j) = (Kmod_CircSTD_mat(i,j) + Kmod_CircSTD_mat(i+1,j) + Kmod_CircSTD_mat(i,j+1) + Kmod_CircSTD_mat(i+1,j+1))/4;
        KmodParamSTD_TMet_mat2(i,j) = (KmodParamSTD_TMet_mat(i,j) + KmodParamSTD_TMet_mat(i+1,j) + KmodParamSTD_TMet_mat(i,j+1) + KmodParamSTD_TMet_mat(i+1,j+1))/4;
        KmodParamSTD_Circ_mat2(i,j) = (KmodParamSTD_Circ_mat(i,j) + KmodParamSTD_Circ_mat(i+1,j) + KmodParamSTD_Circ_mat(i,j+1) + KmodParamSTD_Circ_mat(i+1,j+1))/4;
    end
end

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            Kmod_fwd2(idx, 1) = Kmod_fwd_mat2(i,j);
            Kmod_TMet2(idx, 1) = Kmod_TMet_mat2(i,j);
            Kmod_Circ2(idx, 1) = Kmod_Circ_mat2(i,j);
            KmodParam_TMet2(idx, 1) = KmodParam_TMet_mat2(i,j);
            KmodParam_Circ2(idx, 1) = KmodParam_Circ_mat2(i,j);
        else
            Kmod_fwd2(idx, 1) = Kmod_fwd_mat(i,j);
            Kmod_TMet2(idx, 1) = Kmod_TMet_mat(i,j);
            Kmod_Circ2(idx, 1) = Kmod_Circ_mat(i,j);
            KmodParam_TMet2(idx, 1) = KmodParam_TMet_mat(i,j);
            KmodParam_Circ2(idx, 1) = KmodParam_Circ_mat(i,j);
        end
        idx = idx + 1;
    end
end
idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            KmodSTD_TMet2(idx, 1) = Kmod_TMetSTD_mat2(i,j);
            KmodSTD_Circ2(idx, 1) = Kmod_CircSTD_mat2(i,j);
            KmodParamSTD_TMet2(idx, 1) = KmodParamSTD_TMet_mat2(i,j);
            KmodParamSTD_Circ2(idx, 1) = KmodParamSTD_Circ_mat2(i,j);
        else
            KmodSTD_TMet2(idx, 1) = Kmod_TMetSTD_mat(i,j);
            KmodSTD_Circ2(idx, 1) = Kmod_CircSTD_mat(i,j);
            KmodParamSTD_TMet2(idx, 1) = KmodParamSTD_TMet_mat(i,j);
            KmodParamSTD_Circ2(idx, 1) = KmodParamSTD_Circ_mat(i,j);
        end
        idx = idx + 1;
    end
end

%ModulusTable = [G_vals, K_vals, Kmod_fwd, Kmod_TMet, Kmod_Circ, KmodParam_TMet, KmodParam_Circ, KmodParamSTD_TMet, KmodParamSTD_Circ]
ModulusTable = [G_vals, K_vals, Kmod_fwd2, Kmod_TMet2, Kmod_Circ2, KmodParam_TMet2, KmodParam_Circ2, KmodParamSTD_TMet2, KmodParamSTD_Circ2];

end