function [pl] = plotShapes(DataStruc, subplot_no, RefOrDef)

KG = 'K050_G050';
Strain = DataStruc.(KG).Params_phys.frac;
Wo = DataStruc.(KG).Params_phys.Wo_paper;
Ar = DataStruc.(KG).Params_phys.Ar_paper;
Sigma = DataStruc.(KG).Params_phys.sigma_dimal;
K = DataStruc.(KG).Params_phys.Kmod_dimal;
G = DataStruc.(KG).Params_phys.Gmod_dimal;

%subplot(1,4,subplot_no)
ii = 7;
LW = 2; MS = 7;
Color_Ref = [47, 15, 61]/255;
% Color_Strain95 = [107, 24, 93]/255;
% Color_Strain90 = [168, 40, 96]/255;
% Color_Strain80 = [216, 80, 83]/255;
Color_Strain95 = [138, 29, 99]/255;
Color_Strain90 = [195, 56, 90]/255;
Color_Strain80 = [231, 109, 84]/255;

switch subplot_no
    case 1
        col1 = '(a)';
    case 2
        col1 = '(b)';
    case 3
        col1 = '(c)';
    case 4
        col1 = '(d)';
end

if isequal(RefOrDef, 'Ref')
    pl = plot(DataStruc.(KG).RefSt_fwd.r_dimal, DataStruc.(KG).RefSt_fwd.z_dimal, 'Color', Color_Ref, ...
        'LineWidth', LW); hold on
    xlabel('$r$ [mm]', 'Interpreter', 'latex'); ylabel('$z$ [mm]', 'Interpreter', 'latex');
    title({[strcat(col1, ' Wo=', num2str(Wo), ', Ar=', num2str(Ar), ', $\sigma_{\alpha\beta}$=', num2str(Sigma), ' mN/m, ')],...
        [strcat('$K$=', num2str(K), ' mN/m, $G$=', num2str(G), 'mN/m')]}, 'Interpreter', 'latex');
else
    if Strain == 0.8
        pl = plot(DataStruc.(KG).DefSt_fwd.r_dimal, DataStruc.(KG).DefSt_fwd.z_dimal, 'o:', 'MarkerIndices', 1:ii:length(DataStruc.(KG).DefSt_fwd.r_dimal), ...
            'MarkerEdgeColor', Color_Strain80, 'MarkerFaceColor', Color_Strain80, 'MarkerSize', MS, 'LineWidth', LW, ...
            'Color', Color_Strain80); hold on
    elseif Strain == 0.9
        pl = plot(DataStruc.(KG).DefSt_fwd.r_dimal, DataStruc.(KG).DefSt_fwd.z_dimal, '^:', 'MarkerIndices', 1:ii:length(DataStruc.(KG).DefSt_fwd.r_dimal), ...
            'MarkerEdgeColor', Color_Strain90, 'MarkerFaceColor', Color_Strain90, 'MarkerSize', MS+1, 'LineWidth', LW, ...
            'Color', Color_Strain90); hold on
    elseif Strain == 0.95
        pl = plot(DataStruc.(KG).DefSt_fwd.r_dimal, DataStruc.(KG).DefSt_fwd.z_dimal, 's:', 'MarkerIndices', 1:ii:length(DataStruc.(KG).DefSt_fwd.r_dimal), ...
            'MarkerEdgeColor', Color_Strain95, 'MarkerFaceColor', Color_Strain95, 'MarkerSize', MS+1, 'LineWidth', LW, ...
            'Color', Color_Strain95); hold on
    else
        pl = plot(DataStruc.(KG).DefSt_fwd.r_dimal, DataStruc.(KG).DefSt_fwd.z_dimal, '.', 'k', Color_Strain90, ...
            'LineWidth', LW); hold on
    end
end

set(gca,'DataAspectRatio',[1 1 1], 'XLim', [0, max(DataStruc.(KG).RefSt_fwd.r_dimal) + max(DataStruc.(KG).RefSt_fwd.r_dimal)/10], ...
    'YLim', [min(DataStruc.(KG).RefSt_fwd.z_dimal) + min(DataStruc.(KG).RefSt_fwd.z_dimal)/10, 0], 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', ...
    'TickLabelInterpreter','latex');
