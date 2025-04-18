function [] = getDilParamPlot_v2(DataStruc, VarList, DilParam_table, X, Y, subpl)

plot_lims = [0 2.2];

% make exceptions for:
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.9, stress = Balemans: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Bal90')
    idx = find(DilParam_table(:,1) == 19);
    DilParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.8, stress = Balemans: 
% change G = 49 to G = 50, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal80')
    idx = find(DilParam_table(:,1) == 49);
    DilParam_table(idx,1) = 50;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Pepicelli: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Pep90')
    idx = find(DilParam_table(:,1) == 9);
    DilParam_table(idx,1) = 10;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.9, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Pep90')
    idx = find(DilParam_table(:,1) == 19);
    DilParam_table(idx,1) = 20;
end
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.8, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Pep80')
    idx = find(DilParam_table(:,1) == 19);
    DilParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.95, stress = Balemans: 
% change G = 18 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal95')
    idx = find(DilParam_table(:,1) == 18);
    DilParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Balemans: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Bal90')
    idx = find(DilParam_table(:,1) == 9);
    DilParam_table(idx,1) = 10;
end


Z = NaN(size(X,1));
for k = 1:size(DilParam_table,1)
    if round(DilParam_table(k,2)) == 0 || round(DilParam_table(k,1)) == 0
        continue
    end
    i = find(round(DilParam_table(k,2)) == round(Y(:,1))); % K
    j = find(round(DilParam_table(k,1)) == round(X(1,:))); % G
    Z(i,j) = DilParam_table(k,3);
end
%subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.2 0.1], [0.1 0.1]);
%subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.2 0.2], [0.05 0.1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.1 0.2], [0.05 0.001]);
subplot(1,5,subpl); hold on
MyMap = brewermap(9, 'Reds'); %MyMap = flipud(MyMap);
pl = pcolor(X,Y,Z); pl.EdgeColor = 'none';
%pl.FaceColor = 'interp';
colormap(MyMap);
caxis(plot_lims); %colorbar
ax = gca; ax.TickLength = [0.03, 0.03];
set(gca,'xscale','log','yscale','log', 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', 'layer', 'top', 'TickLabelInterpreter','latex');
grid on
ax.GridColor = [0 0 0]; ax.GridAlpha = 1; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 1; ax.MinorGridLineStyle = '--';
xlabel('Shear modulus, G [mN/m]', 'interpreter', 'latex');
if subpl == 1
    ylabel('Dilatational modulus, K [mN/m]', 'interpreter', 'latex');
else
    set(gca, 'YTickLabel',[]);
end

xl = xlim; yl = ylim;
plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5);  % Left & Lower Axes
plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5);  % Right & Upper Axes
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

%TitleLabel = matlab.lang.makeValidName(strcat('DAP: ', DilParam));
ttl = title(strcat(col, 'Wo = ', num2str(DataStruc.(VarList{i,j}).Params_phys.Wo_paper), ...
    ', Ar = ', num2str(DataStruc.(VarList{i,j}).Params_phys.Ar_paper)), 'interpreter', 'latex');
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
hold on

if subpl == 4
    subplot(1,4,5)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    hCB.Position = [0.76, 0.2, 0.015, 0.61];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
end

Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
StrainMeasure = DataStruc.(VarList{1,1}).Params_phys.strainmeasure;
Strain = round(100*(1 - DataStruc.(VarList{1,1}).Params_phys.frac));
sgtitle(strcat('Dilatational Anisotropy Parameter (DAP): $\Delta A / A$=', num2str(Strain), ...
    '\%, $\sigma_0$=', num2str(Sigma),  ' mN/m, Strain Measure=', StrainMeasure, '\hspace{3.5cm}'), ...
    'interpreter', 'latex', 'FontSize', 18);

end