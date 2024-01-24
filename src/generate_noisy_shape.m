function [rr_noise,zz_noise] = generate_noisy_shape(vars_sol,vars_num,Nsample,sigma_noise)

    Nsample_full = 2*Nsample-1;

    [s_plot,r_plot,z_plot,normals_plot(:,1),normals_plot(:,2)] = ...
        interpolate_on_uniform_grid(vars_num, ...
        [vars_sol.r, vars_sol.z, vars_sol.normals], Nsample); 
    
    [s_plot_full,r_plot_full,z_plot_full,normals_plot_full] = ...
        mirror_shape(s_plot,r_plot,z_plot,normals_plot);
    
    tmp = normrnd(0,sigma_noise,[Nsample_full,1]);
    rr_noise = zeros(Nsample_full,1); zz_noise = rr_noise;
    for i=1:Nsample_full
        rr_noise(i) = r_plot_full(i) + tmp(i)*normals_plot_full(i,1);
        zz_noise(i) = z_plot_full(i) + tmp(i)*normals_plot_full(i,2);
    end

end

