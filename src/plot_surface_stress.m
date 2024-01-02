function [] =     plot_surface_stress(s, sigmas, sigmap, fig_num)

    % plot the surface stresses
    figure(fig_num); hold on
    h(1) = plot(s,sigmas,'LineWidth',2); hold on
    h(2) = plot(s,sigmap,'LineWidth',2);
    xlabel('s','FontSize',32);
    ylabel('\sigma','FontSize',32);
    legend(h,'\sigma_s','\sigma_\phi','FontSize',24,'Location','northwest');
    xlim([0,s(end)])
    ax = gca; ax.FontSize = 24;

end

