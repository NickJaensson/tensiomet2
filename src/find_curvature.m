function [kappas,kappap] = find_curvature(vars_sol,vars_num,make_plot)

    % determine the curvatures
    % NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
    % solved here by taking kappap(0) = kappas(0)
    kappas = vars_num.Ds*vars_sol.psi;
    kappap = kappas;
    kappap(2:end) = sin(vars_sol.psi(2:end))./vars_sol.r(2:end);
    
    if make_plot
        % plot the curvatures versus the z-coordinate
        figure;
        plot(vars_sol.z,kappas,'LineWidth',2); hold on
        plot(vars_sol.z,kappap,'LineWidth',2); hold on
        plot(vars_sol.z,kappap+kappas,'LineWidth',2);
        xlabel('z','FontSize',32);
        ylabel('\kappa','FontSize',32);
        legend('\kappa_s','\kappa_\phi','\kappa_s+\kappa_\phi', ...
            'FontSize',24,'Location','northwest');
        xlim([vars_sol.z(1),0])
        ax = gca; ax.FontSize = 24;
    end

end

