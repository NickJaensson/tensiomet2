function [] = compareStrainStressProfilesPlot_NeoHookean(InputDataStruc, InputVarList, subpl)

% All elements in InputDataStruc must have the same strain, Wo, and Ar values

% Get reference values from 1st entry in DataStruc
strain = InputDataStruc.(InputVarList{1,1}).Params_phys.frac;
sigma0 = InputDataStruc.(InputVarList{1,1}).Params_phys.sigma_dimal;
StrainMeasure = InputDataStruc.(InputVarList{1,1}).Params_phys.strainmeasure;
strainlabel = num2str(round((1 - InputDataStruc.(InputVarList{1,1}).Params_phys.frac)*100));
colors = [0, 204, 0; 0, 102, 255; 255, 51, 153]/255; % green, blue, pink
markers = ['s', '^', 'o'];
marker_idx = [1:10:length(InputDataStruc.(InputVarList{1,1}).DefSt_fwd.s)];
LW = 2; MS = 7;

% Change axes labels, depending on the strain value
switch strain
    case 0.95
        ylim_lt = [0.896 1.005];
        if sigma0 == 60
            ylim_lb = [55 60.2];
        elseif sigma0 == 20
            ylim_lb = [16 20.3];
        end
    case 0.90
        ylim_lt = [0.795 1.01];
        if sigma0 == 60
            ylim_lb = [50 60.4];
        elseif sigma0 == 20
            ylim_lb = [11 20.5];
        end
    case 0.80
        ylim_lt = [0.60 1.02];
        if sigma0 == 60
            ylim_lb = [35 60.8];
        elseif sigma0 == 20
            ylim_lb = [0 21];
        end
end

switch subpl
    case 1
        col = '(a) \,';
    case 2
        col = '(b) \,';
    case 3
        col = '(c) \,';
    case 4
        col = '(d) \,';
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.02], [0.075 0.075], [0.07 0.01]);
subplot(2,4,subpl); hold on

% 1. Plot strains
for i = 1:size(InputVarList,2)
    K_vals(i, 1) = InputDataStruc.(InputVarList{1,i}).Params_phys.Kmod_dimal;
    G_vals(i, 1) = InputDataStruc.(InputVarList{1,i}).Params_phys.Gmod_dimal;
    % normalize dimensionless arc length, such that s ranges from 0 to 1
    s_norm = InputDataStruc.(InputVarList{1,i}).DefSt_fwd.s./InputDataStruc.(InputVarList{1,i}).DefSt_fwd.s(end);
    % find dilatational strain along arc length
    J = InputDataStruc.(InputVarList{1,i}).DefSt_fwd.lams.*InputDataStruc.(InputVarList{1,i}).DefSt_fwd.lamp;
    strain_dil = log(J)./J;
    % find shear strain along arc length
    lams2Inv = 1./((InputDataStruc.(InputVarList{1,i}).DefSt_fwd.lams).^2);
    lamp2Inv = 1./((InputDataStruc.(InputVarList{1,i}).DefSt_fwd.lamp).^2);
    strain_shr = 0.5*(lamp2Inv - lams2Inv);
    % make plots
    plot(s_norm, strain_dil, '-',  'LineWidth', LW, 'Color', colors(i,:), 'Marker', markers(i), 'MarkerIndices', marker_idx, 'MarkerFaceColor', colors(i,:), 'MarkerSize', MS);
    plot(s_norm, strain_shr, '--', 'LineWidth', LW, 'Color', colors(i,:), 'Marker', markers(i), 'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'MarkerSize', MS);
    p(i) = plot(s_norm(1), strain_dil(1), 'LineStyle', 'none', 'Marker', markers(i), 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', MS);
end

Wo = InputDataStruc.(InputVarList{1,1}).Params_phys.Wo_paper;
Ar = InputDataStruc.(InputVarList{1,1}).Params_phys.Ar_paper;
set(gca, 'XTickLabel',[]);
%set(gca,'ylim', ylim_lt, 'LineWidth', 1, 'FontSize', 14, 'TickDir','out', 'TickLabelInterpreter','latex');
set(gca,'LineWidth', 1, 'FontSize', 14, 'TickDir','out', 'TickLabelInterpreter','latex');
grid on; ax = gca;
ax.GridColor = [0 0 0]; ax.GridAlpha = 0.4; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 0.4; ax.MinorGridLineStyle = '--';
ttl = title(strcat({strcat(col, 'Wo=', num2str(Wo), ', Ar=', num2str(Ar), ', $\sigma_0$=', num2str(sigma0), ' mN/m,') ; ...
    strcat('\qquad $\Delta A/A_0$=', strainlabel, '\%, Strain measure=', StrainMeasure)}), 'interpreter', 'latex');
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
if subpl == 1
    ylabel({'Dilatational (\textemdash) and shear (- -)' ; 'deformation fractions, [ ]'}, 'interpreter', 'latex');
    %ylabel({'$\lambda_s\lambda_p$ (\textemdash) and $\lambda_s/\lambda_p$ (- -), [ ]'}, 'interpreter', 'latex');
else
    set(gca, 'YTickLabel',[]);
end

% 2. Plot stresses
subplot(2,4,subpl+4); hold on
for i = 1:size(InputVarList,2)
    % normalize dimensionless arc length, such that s ranges from 0 to 1
    s_norm = InputDataStruc.(InputVarList{1,i}).DefSt_fwd.s./InputDataStruc.(InputVarList{1,i}).DefSt_fwd.s(end);
    % find meridional stress along arc length
    stress_s = InputDataStruc.(InputVarList{1,i}).DefSt_fwd.sigmas_dimal;
    % find azimuthal strain along arc length
    stress_p = InputDataStruc.(InputVarList{1,i}).DefSt_fwd.sigmap_dimal;
    % make plots
    plot(s_norm, stress_s,  '-', 'LineWidth', LW, 'Color', colors(i,:), 'Marker', markers(i), 'MarkerIndices', marker_idx, 'MarkerFaceColor', colors(i,:), 'MarkerSize', MS);
    plot(s_norm, stress_p, '--', 'LineWidth', LW, 'Color', colors(i,:), 'Marker', markers(i), 'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'MarkerSize', MS);
    p(i) = plot(s_norm(1), stress_s(1), 'LineStyle', 'none', 'Marker', markers(i), 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', MS+2);
end
xlabel('Dimensionless arc length, s [ ]', 'interpreter', 'latex');

if subpl == 1
    ylabel({'Meridional (\textemdash) and azimuthal (- -)' ; 'surface stresses [mN/m]'}, 'interpreter', 'latex');
    %ylabel({'$\sigma_s$ (\textemdash) and $\sigma_\phi$ (- -), [mN/m]'}, 'interpreter', 'latex');
elseif subpl == 4
    set(gca, 'YTickLabel',[]);
    labels1 = {strcat('\qquad', num2str(K_vals(1)), '\, \qquad'), strcat('\qquad', num2str(K_vals(2)), '\qquad'), strcat('\qquad', num2str(K_vals(3)), '\qquad')};
    labels2 = {strcat('\qquad', num2str(G_vals(1)), '\quad'), strcat('\qquad', num2str(G_vals(2)), '\quad'), strcat('\qquad', num2str(G_vals(3)), '\quad')};
    labels = cellfun(@(x,y) {sprintf('%1s%1s', x, y)}, labels1, labels2);
    l = legend(p, labels, 'interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
    titles = {'\quad \qquad', '$K$ (mN/m) \,', '$G$ (mN/m)'};
    titles = sprintf('%1s%1s%1s', titles{:});
    l.Title.String = titles;
else
    set(gca, 'YTickLabel',[]);
end
set(gca,'ylim', ylim_lb, 'LineWidth', 1, 'FontSize', 14, 'TickDir', 'out', 'TickLabelInterpreter','latex');
grid on; ax = gca;
ax.GridColor = [0 0 0]; ax.GridAlpha = 0.4; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 0.4; ax.MinorGridLineStyle = '--';


end
