function [rr_noise, zz_noise] = ...
    generate_noisy_shape(vars_sol, vars_num, Nsample, sigma_noise)
    % GENERATE_NOISY_SHAPE Generates full shape w/ added noise from forward 
    % solution.
    %
    % INPUTS:
    %   vars_sol     - Structure with solution variables (r, z, normals)
    %   vars_num     - Structure with numerical parameters
    %   Nsample      - Number of points of forward solution ("half" shape)
    %   sigma_noise  - Standard deviation of the noise to add
    %
    % OUTPUTS:
    %   rr_noise     - Noisy radial coordinates
    %   zz_noise     - Noisy axial coordinates

    Nsample_full = 2*Nsample-1;

    [s_plot,r_plot,z_plot,normals_plot(:,1),normals_plot(:,2)] = ...
        interpolate_on_uniform_grid(vars_num, ...
        [vars_sol.r, vars_sol.z, vars_sol.normals], Nsample); 
    
    [s_plot_full,r_plot_full,z_plot_full,normals_plot_full] = ...
        mirror_shape(s_plot,r_plot,z_plot,normals_plot);
    
    % add noise on all points, except the first and last point
    tmp = normrnd(0,sigma_noise,[Nsample_full,1]);
    rr_noise = zeros(Nsample_full,1); zz_noise = rr_noise;
    for i=2:Nsample_full-1
        rr_noise(i) = r_plot_full(i) + tmp(i)*normals_plot_full(i,1);
        zz_noise(i) = z_plot_full(i) + tmp(i)*normals_plot_full(i,2);
    end
    rr_noise(1) = r_plot_full(1);
    zz_noise(1) = z_plot_full(1);
    rr_noise(Nsample_full) = r_plot_full(Nsample_full);
    zz_noise(Nsample_full) = z_plot_full(Nsample_full);

end

