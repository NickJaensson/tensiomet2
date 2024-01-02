function [normals] = get_normals(vars_sol, vars_num)

    % determine the normal vectors
    normals(:,1) = vars_num.Ds*vars_sol.z; % r-component of the normals
    normals(:,2) = -vars_num.Ds*vars_sol.r; % z-component of the normals
    for i=1:size(normals,2)
        normals(i,:) = normals(i,:)/norm(normals(i,:));
    end

end

