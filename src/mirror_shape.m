function [s_full,r_full,z_full,normals_full] = mirror_shape(s,r,z,normals)

    tmp = s-s(end);
    s_full = [tmp(1:end-1)',s']';
    
    r_full = [-flipud(r)',r(2:end)']';
    
    z_full = [flipud(z)',z(2:end)']';
    
    normals_full(:,1) = [-flipud(normals(:,1))',normals(2:end,1)']';
    normals_full(:,2) = [flipud(normals(:,2))',normals(2:end,2)']';

end

