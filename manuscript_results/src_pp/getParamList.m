function [VarList] = getParamList(K_vec, G_vec, K_str, G_str)

for i = 1:length(K_vec)
    for j = 1:length(G_vec)
        VarList{i,j} = matlab.lang.makeValidName(strcat('K', K_str{i}, '_G', G_str{j}));
    end
end

end