function [] =     plot_surface_stress(vars_num,vars_sol)

    % plot the surface stresses
    figure;
    plot(vars_num.s,vars_sol.sigmas,'LineWidth',2); hold on
    plot(vars_num.s,vars_sol.sigmap,'LineWidth',2);
    xlabel('s','FontSize',32);
    ylabel('\sigma','FontSize',32);
    legend('\sigma_s','\sigma_\phi','FontSize',24,'Location','northwest');
    xlim([0,vars_num.s(end)])
    ax = gca; ax.FontSize = 24;

end

