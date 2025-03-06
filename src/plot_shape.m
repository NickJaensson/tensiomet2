function [] = plot_shape(r, z, fig_num)
    % PLOT_SHAPE Plots the shape profile in the (r, z) plane.
    %
    % INPUTS:
    %   r       - Radial coordinates.
    %   z       - Axial coordinates.
    %   fig_num - Figure number for plotting.

    figure(fig_num); hold on
    plot(r, z,'-o');
    set(gca,'DataAspectRatio',[1 1 1])

end
