function [] =     plot_curvature(z, kappas, kappap, fig_num)

    % plot the curvatures versus the z-coordinate
    figure(fig_num); hold on
    plot(z,kappas,'LineWidth',2); hold on
    plot(z,kappap,'LineWidth',2); hold on
    plot(z,kappap+kappas,'LineWidth',2);
    xlabel('z','FontSize',32);
    ylabel('\kappa','FontSize',32);
    legend('\kappa_s','\kappa_\phi','\kappa_s+\kappa_\phi', ...
        'FontSize',24,'Location','northwest');
    xlim([z(1),0])
    ax = gca; ax.FontSize = 24;

end

