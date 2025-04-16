function [] =     plot_surface_deformation(s, lams, lamp, fig_num, style)
    % PLOT_SURFACE_DEFORMATION Plots the surface deformations:
    %     \lambda_s * \lambda_\phi (red curve)
    %     \lambda_s / \lambda_\phi (blue curve)
    %
    % INPUTS:
    %   s       - Arclength coordinates.
    %   lams    - Meridional stretch ratio (\lambda_s).
    %   lamp    - Circumferential stretch ratio (\lambda_\phi).
    %   fig_num - Figure number for plotting.
    %   style   - Optional: Line style for plotting (string, e.g., '--').

    if nargin < 5
        style = '-';
    end

    % plot the surface deformations
    figure(fig_num); hold on
    h(1) = plot(s,lams.*lamp,'LineWidth',2,'Color','r','LineStyle',style); 
    hold on
    h(2) = plot(s,lams./lamp,'LineWidth',2,'Color','b','LineStyle',style);
    xlabel('s [-]','FontSize',32);
    ylabel('Deformation, [-]','FontSize',32);
    legend(h,'\lambda_s * \lambda_\phi','\lambda_s / \lambda_\phi', ...
        'FontSize',24,'Location','northwest');
    xlim([0,s(end)])
    ax = gca; ax.FontSize = 24;

end
