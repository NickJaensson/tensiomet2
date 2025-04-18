function [ApexStrainTable] = getStrainApexNonhomogen(DataStruc, VarList)

ApexStrain_fwd  = zeros(numel(VarList), 1);
ApexStrain_TMet = zeros(numel(VarList), 1);
ApexStrain_Circ = zeros(numel(VarList), 1);
StrainParam_TMet = zeros(numel(VarList), 1);
StrainParam_Circ = zeros(numel(VarList), 1);
K_vals = zeros(numel(VarList), 1);
G_vals = zeros(numel(VarList), 1);

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        K_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        G_vals(idx, 1) = DataStruc.(VarList{i,j}).Params_phys.Gmod_dimal;
        % Get dilatational strain at the apex (lams*lamp)
        ApexStrain_fwd(idx,1)  = DataStruc.(VarList{i,j}).DefSt_fwd.lams(1,1)*DataStruc.(VarList{i,j}).DefSt_fwd.lamp(1,1);
        ApexStrain_TMet(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.lams(1,1)*DataStruc.(VarList{i,j}).DefSt_inv.lamp(1,1);
        ApexStrain_Circ(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.area_dimal/DataStruc.(VarList{i,j}).RefSt_inv.area_dimal;
        ApexStrain_fwd_mat(i,j) = ApexStrain_fwd(idx,1);
        ApexStrain_TMet_mat(i,j) = ApexStrain_TMet(idx,1);
        ApexStrain_Circ_mat(i,j) = ApexStrain_Circ(idx,1);
        % strain param = ((1-StrainModel) - (1-StrainReal)) / (1-StrainReal)
        StrainParam_TMet(idx,1) = ((1-ApexStrain_TMet(idx,1)) - (1-ApexStrain_fwd(idx,1)))/(1-ApexStrain_fwd(idx,1));
        StrainParam_Circ(idx,1) = ((1-ApexStrain_Circ(idx,1)) - (1-ApexStrain_fwd(idx,1)))/(1-ApexStrain_fwd(idx,1));
        StrainParam_TMet_mat(i,j) = StrainParam_TMet(idx,1);
        StrainParam_Circ_mat(i,j) = StrainParam_Circ(idx,1);
        idx = idx + 1;
    end
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
frac = DataStruc.(VarList{1,1}).Params_phys.frac;
for i = 1:size(VarList,1)-1
    for j = 1:size(VarList, 2)-1
        if  Sigma == 20 && frac == 0.8  %excption for inaccessible geometries, average less cells
            if i == size(VarList,1)-1 
                ApexStrain_fwd_mat2(i,j) = ApexStrain_fwd_mat(i,j);
                ApexStrain_TMet_mat2(i,j) = ApexStrain_TMet_mat(i,j);
                ApexStrain_Circ_mat2(i,j) = ApexStrain_Circ_mat(i,j);
                StrainParam_TMet_mat2(i,j) = StrainParam_TMet_mat(i,j);
                StrainParam_Circ_mat2(i,j) = StrainParam_Circ_mat(i,j);
            else
                ApexStrain_fwd_mat2(i,j) = (ApexStrain_fwd_mat(i,j) + ApexStrain_fwd_mat(i+1,j) + ApexStrain_fwd_mat(i,j+1) + ApexStrain_fwd_mat(i+1,j+1))/4;
                ApexStrain_TMet_mat2(i,j) = (ApexStrain_TMet_mat(i,j) + ApexStrain_TMet_mat(i+1,j) + ApexStrain_TMet_mat(i,j+1) + ApexStrain_TMet_mat(i+1,j+1))/4;
                ApexStrain_Circ_mat2(i,j) = (ApexStrain_Circ_mat(i,j) + ApexStrain_Circ_mat(i+1,j) + ApexStrain_Circ_mat(i,j+1) + ApexStrain_Circ_mat(i+1,j+1))/4;
                StrainParam_TMet_mat2(i,j) = (StrainParam_TMet_mat(i,j) + StrainParam_TMet_mat(i+1,j) + StrainParam_TMet_mat(i,j+1) + StrainParam_TMet_mat(i+1,j+1))/4;
                StrainParam_Circ_mat2(i,j) = (StrainParam_Circ_mat(i,j) + StrainParam_Circ_mat(i+1,j) + StrainParam_Circ_mat(i,j+1) + StrainParam_Circ_mat(i+1,j+1))/4;
            end
        else
            ApexStrain_fwd_mat2(i,j) = (ApexStrain_fwd_mat(i,j) + ApexStrain_fwd_mat(i+1,j) + ApexStrain_fwd_mat(i,j+1) + ApexStrain_fwd_mat(i+1,j+1))/4;
            ApexStrain_TMet_mat2(i,j) = (ApexStrain_TMet_mat(i,j) + ApexStrain_TMet_mat(i+1,j) + ApexStrain_TMet_mat(i,j+1) + ApexStrain_TMet_mat(i+1,j+1))/4;
            ApexStrain_Circ_mat2(i,j) = (ApexStrain_Circ_mat(i,j) + ApexStrain_Circ_mat(i+1,j) + ApexStrain_Circ_mat(i,j+1) + ApexStrain_Circ_mat(i+1,j+1))/4;
            StrainParam_TMet_mat2(i,j) = (StrainParam_TMet_mat(i,j) + StrainParam_TMet_mat(i+1,j) + StrainParam_TMet_mat(i,j+1) + StrainParam_TMet_mat(i+1,j+1))/4;
            StrainParam_Circ_mat2(i,j) = (StrainParam_Circ_mat(i,j) + StrainParam_Circ_mat(i+1,j) + StrainParam_Circ_mat(i,j+1) + StrainParam_Circ_mat(i+1,j+1))/4;
        end
    end
end

idx = 1;
for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        if i < size(VarList,1) && j < size(VarList,2)
            ApexStrain_fwd2(idx, 1) = ApexStrain_fwd_mat2(i,j);
            ApexStrain_TMet2(idx, 1) = ApexStrain_TMet_mat2(i,j);
            ApexStrain_Circ2(idx, 1) = ApexStrain_Circ_mat2(i,j);
            StrainParam_TMet2(idx, 1) = StrainParam_TMet_mat2(i,j);
            StrainParam_Circ2(idx, 1) = StrainParam_Circ_mat2(i,j);
        else
            ApexStrain_fwd2(idx, 1) = ApexStrain_fwd_mat(i,j);
            ApexStrain_TMet2(idx, 1) = ApexStrain_TMet_mat(i,j);
            ApexStrain_Circ2(idx, 1) = ApexStrain_Circ_mat(i,j);
            StrainParam_TMet2(idx, 1) = StrainParam_TMet_mat(i,j);
            StrainParam_Circ2(idx, 1) = StrainParam_Circ_mat(i,j);
        end
        idx = idx + 1;
    end
end

%ApexStrainTable = [G_vals, K_vals, ApexStrain_fwd, ApexStrain_TMet, ApexStrain_Circ, StrainParam_TMet, StrainParam_Circ]
ApexStrainTable = [G_vals, K_vals, ApexStrain_fwd2, ApexStrain_TMet2, ApexStrain_Circ2, StrainParam_TMet2, StrainParam_Circ2];

end