function [] = plot_shape(r, z, fig_num)

    figure(fig_num); hold on
    plot(r, z,'-o');
    set(gca,'DataAspectRatio',[1 1 1])

end

