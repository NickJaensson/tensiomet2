function [DataStruc] = assembleResults(Wo, Ar, Sigma, StrainMeasure, Strain, K_vec, G_vec, VarList)

for i = 1:length(K_vec)
    for j = 1:length(G_vec)
        FolderName = strcat('Results/Wo=', num2str(Wo), '_Ar=', num2str(Ar), ...
            '_sigma0=', num2str(Sigma), '_K=', num2str(K_vec(i)), '_G=', num2str(G_vec(j)), ...
            '_StrainM=', StrainMeasure, '_Strain=', num2str(Strain));
        DataStruc.(VarList{i,j}) = loadResults(FolderName);
    end
end

end