function [] = plot_shape(z_plot,r_plot,vars_sol)

    % plot the shape of the drop on the plotting grid

    % if ~isfield( vars_sol, 'r_star' )
         figure; hold on
    % end
    scatter(r_plot',z_plot','b');
    plot(r_plot',z_plot','b');

    % rescale the plot
    % if isfield( vars_sol, 'r_star' )
    %     rmax = max([vars_sol.r_star',r_plot']);
    %     zmin = min([vars_sol.z_star',z_plot']);
    %     xlim([0 1.2*rmax]);
    %     ylim([1.2*zmin 0]);
    % end
    set(gca,'DataAspectRatio',[1 1 1])

end

