function [nnormals] = get_normals(vars_sol, vars_num, s_plot)
    % determine the normal vectors
    normals(:,1) = vars_num.Ds*vars_sol.z; % r-component of the normals
    normals(:,2) = -vars_num.Ds*vars_sol.r; % z-component of the normals
    for i=1:size(normals,2)
        normals(i,:) = normals(i,:)/norm(normals(i,:));
    end
    
    nnormals(:,1) = interp1(vars_num.s,normals(:,1),s_plot,'pchip');
    nnormals(:,2) = interp1(vars_num.s,normals(:,2),s_plot,'pchip');
    for i=1:size(nnormals,2)
        nnormals(i,:) = nnormals(i,:)/norm(nnormals(i,:));
    end
end

