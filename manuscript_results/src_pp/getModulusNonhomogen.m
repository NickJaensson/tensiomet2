function [ModulusTable] = getModulusNonhomogen(DataStruc, VarList)

Kmod_fwd  = zeros(numel(VarList), 1);
Kmod_TMet = zeros(numel(VarList), 1);
Kmod_Circ = zeros(numel(VarList), 1);
KmodParam_TMet = zeros(numel(VarList), 1);
KmodParam_Circ = zeros(numel(VarList), 1);
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
        Kmod_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.Kmod_dimal;
        Kmod_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.KmodCircleFit_dimal;
        Kmod_fwd_mat(i,j) = Kmod_fwd(idx,1);
        Kmod_TMet_mat(i,j) = Kmod_TMet(idx,1);
        Kmod_Circ_mat(i,j) = Kmod_Circ(idx,1);
        % modulus param = (ModulusModel - ModulusReal) / ModulusReal
        KmodParam_TMet(idx,1) = (Kmod_TMet(idx,1) - Kmod_fwd(idx,1))/Kmod_fwd(idx,1);
        KmodParam_Circ(idx,1) = (Kmod_Circ(idx,1) - Kmod_fwd(idx,1))/Kmod_fwd(idx,1);
        KmodParam_TMet_mat(i,j) = KmodParam_TMet(idx,1);
        KmodParam_Circ_mat(i,j) = KmodParam_Circ(idx,1);
        idx = idx + 1;
    end
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
frac = DataStruc.(VarList{1,1}).Params_phys.frac;
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        if  Sigma == 20 && frac == 0.8  %excption for inaccessible geometries, average less cells
            if i == size(VarList,1)-1 
                Kmod_fwd_mat2(i,j) = Kmod_fwd_mat(i,j);
                Kmod_TMet_mat2(i,j) = Kmod_TMet_mat(i,j);
                Kmod_Circ_mat2(i,j) = Kmod_Circ_mat(i,j);
                KmodParam_TMet_mat2(i,j) = KmodParam_TMet_mat(i,j);
                KmodParam_Circ_mat2(i,j) = KmodParam_Circ_mat(i,j);
            else
                Kmod_fwd_mat2(i,j) = (Kmod_fwd_mat(i,j) + Kmod_fwd_mat(i+1,j) + Kmod_fwd_mat(i,j+1) + Kmod_fwd_mat(i+1,j+1))/4;
                Kmod_TMet_mat2(i,j) = (Kmod_TMet_mat(i,j) + Kmod_TMet_mat(i+1,j) + Kmod_TMet_mat(i,j+1) + Kmod_TMet_mat(i+1,j+1))/4;
                Kmod_Circ_mat2(i,j) = (Kmod_Circ_mat(i,j) + Kmod_Circ_mat(i+1,j) + Kmod_Circ_mat(i,j+1) + Kmod_Circ_mat(i+1,j+1))/4;
                KmodParam_TMet_mat2(i,j) = (KmodParam_TMet_mat(i,j) + KmodParam_TMet_mat(i+1,j) + KmodParam_TMet_mat(i,j+1) + KmodParam_TMet_mat(i+1,j+1))/4;
                KmodParam_Circ_mat2(i,j) = (KmodParam_Circ_mat(i,j) + KmodParam_Circ_mat(i+1,j) + KmodParam_Circ_mat(i,j+1) + KmodParam_Circ_mat(i+1,j+1))/4;
            end
        else
            Kmod_fwd_mat2(i,j) = (Kmod_fwd_mat(i,j) + Kmod_fwd_mat(i+1,j) + Kmod_fwd_mat(i,j+1) + Kmod_fwd_mat(i+1,j+1))/4;
            Kmod_TMet_mat2(i,j) = (Kmod_TMet_mat(i,j) + Kmod_TMet_mat(i+1,j) + Kmod_TMet_mat(i,j+1) + Kmod_TMet_mat(i+1,j+1))/4;
            Kmod_Circ_mat2(i,j) = (Kmod_Circ_mat(i,j) + Kmod_Circ_mat(i+1,j) + Kmod_Circ_mat(i,j+1) + Kmod_Circ_mat(i+1,j+1))/4;
            KmodParam_TMet_mat2(i,j) = (KmodParam_TMet_mat(i,j) + KmodParam_TMet_mat(i+1,j) + KmodParam_TMet_mat(i,j+1) + KmodParam_TMet_mat(i+1,j+1))/4;
            KmodParam_Circ_mat2(i,j) = (KmodParam_Circ_mat(i,j) + KmodParam_Circ_mat(i+1,j) + KmodParam_Circ_mat(i,j+1) + KmodParam_Circ_mat(i+1,j+1))/4;
        end
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

%ModulusTable = [G_vals, K_vals, Kmod_fwd, Kmod_TMet, Kmod_Circ, KmodParam_TMet, KmodParam_Circ];
ModulusTable = [G_vals, K_vals, Kmod_fwd2, Kmod_TMet2, Kmod_Circ2, KmodParam_TMet2, KmodParam_Circ2];

end