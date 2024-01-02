function [] = plot_shape(vars_sol, fig_num)

    % plot the shape of the drop on the plotting grid

    figure(fig_num); hold on
    plot(vars_sol.r,vars_sol.z,'-o');

    % rescale the plot
    % if isfield( vars_sol, 'r_star' )
    %     rmax = max([vars_sol.r_star',r_plot']);
    %     zmin = min([vars_sol.z_star',z_plot']);
    %     xlim([0 1.2*rmax]);
    %     ylim([1.2*zmin 0]);
    % end
    set(gca,'DataAspectRatio',[1 1 1])

end

