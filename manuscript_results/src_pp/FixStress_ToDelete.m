function [DataStruc] = ComputeRelevantQuantities_Noise(DataStruc, VarList, NumSimuls)

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        sigmas_apex_vec = zeros(NumSimuls,1);
        for kk = 1:NumSimuls
            sigmas_apex_vec(kk) = DataStruc.(VarList{i,j}).DefSt_inv{kk}.sigmas(1);
        end
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Avg = mean(sigmas_apex_vec);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Std = std(sigmas_apex_vec);
        Redim = DataStruc.(VarList{i,j}).Params_phys.deltarho_dimal*DataStruc.(VarList{i,j}).Params_phys.grav_dimal*(DataStruc.(VarList{i,j}).Params_phys.rneedle_dimal)^2;
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Avg_dimal = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Avg*Redim;
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Std_dimal = DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.sigmas_apex_Std*Redim;
    end
end

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        area_dimal_ref_vec = zeros(NumSimuls,1);
        for kk = 1:NumSimuls
            area_dimal_ref_vec(kk) = DataStruc.(VarList{i,j}).RefSt_inv{kk}.area_dimal;
        end
        DataStruc.(VarList{i,j}).RefSt_inv{NumSimuls+1}.area_Avg_dimal = mean(area_dimal_ref_vec);
        DataStruc.(VarList{i,j}).RefSt_inv{NumSimuls+1}.area_Std_dimal = std(area_dimal_ref_vec);
    end
end

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        apex_strain_Circ_vec = zeros(NumSimuls,1);
        for kk = 1:NumSimuls
            apex_strain_Circ_vec(kk) = (DataStruc.(VarList{i,j}).DefSt_inv{kk}.area_dimal)/DataStruc.(VarList{i,j}).RefSt_inv{kk}.area_dimal;
        end
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Apex_DilStrain_Circ_Avg = mean(apex_strain_Circ_vec);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Apex_DilStrain_Circ_Std = std(apex_strain_Circ_vec);
    end
end

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        apex_strain_vec = zeros(NumSimuls,1);
        for kk = 1:NumSimuls
            apex_strain_vec(kk) = DataStruc.(VarList{i,j}).DefSt_inv{kk}.lamp(1,1)*DataStruc.(VarList{i,j}).DefSt_inv{kk}.lams(1,1);
        end
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Apex_DilStrain_Avg = mean(apex_strain_vec);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.Apex_DilStrain_Std = std(apex_strain_vec);
    end
end

end

