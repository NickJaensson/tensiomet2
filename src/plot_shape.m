function [] = plot_shape(z_plot,r_plot)

    % plot the shape of the drop on the plotting grid
    figure; hold on
    scatter(r_plot',z_plot','b');
    plot(r_plot',z_plot','b');
    set(gca,'DataAspectRatio',[1 1 1])

end

