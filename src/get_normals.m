function [normals] = get_normals(vars_sol, vars_num)
    % GET_NORMALS Computes the unit normal vectors for the shape.
    %
    % INPUTS:
    %   vars_sol - Structure with solution variables
    %   vars_num - Structure with numerical parameters
    %
    % OUTPUTS:
    %   normals  - Normalized normal vectors (r- and z-components)

    % determine the normal vectors
    normals(:,1) = vars_num.Ds*vars_sol.z; % r-component of the normals
    normals(:,2) = -vars_num.Ds*vars_sol.r; % z-component of the normals
    for i=1:size(normals,1)
        normals(i,:) = normals(i,:)/norm(normals(i,:));
    end

end
