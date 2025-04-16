function [DataStruc] = assembleResults_Noise(Wo, Ar, Sigma, StrainMeasure, Strain, K_vec, G_vec, DimlessNoise, NumCheb, NumSimuls, VarList)

for i = 1:length(K_vec)
    for j = 1:length(G_vec)
        FolderName = strcat('Results_Noise/Wo=', num2str(Wo), '_Ar=', num2str(Ar), ...
            '_sigma0=', num2str(Sigma), '_K=', num2str(K_vec(i)), '_G=', num2str(G_vec(j)), ...
            '_StrainM=', StrainMeasure, '_Strain=', num2str(Strain), '_Noise=', num2str(DimlessNoise), ...
            '_NumChebPolys=', num2str(NumCheb));
        DataStruc.(VarList{i,j}) = loadResults_Noise(FolderName, NumSimuls);
    end
end

end