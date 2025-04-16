function [] = getModulusNonhomogenPlotSTD_Noise(DataStruc, VarList, ModulusTable, X, Y, subpl)

% ModulusTable: | G | K | Kmod_fwd | Kmod_TMet | Kmod_Circ | KmodParam_TMet | KmodParam_Circ | KmodParam_TMet_STD | KmodParam_Circ_STD | 

plot_lims = [0 1];
strainlabel = num2str(round((1 - DataStruc.(VarList{1,1}).Params_phys.frac)*100));
Wo = DataStruc.(VarList{1,1}).Params_phys.Wo_paper;
Ar = DataStruc.(VarList{1,1}).Params_phys.Ar_paper;
Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
Strain = round(100*(1 - DataStruc.(VarList{1,1}).Params_phys.frac));
switch subpl
    case 1
        col1 = '(a) \,';
        col2 = '(e) \,';
    case 2
        col1 = '(b) \,';
        col2 = '(f) \,';
    case 3
        col1 = '(c) \,';
        col2 = '(g) \,';
    case 4
        col1 = '(d) \,';
        col2 = '(h) \,';
    otherwise
        col1 = ' ';
        col2 = ' ';
end
TMet_subpltitle = strcat(col1, 'SFE: Wo=', num2str(Wo), ', Ar=', num2str(Ar));
Circ_subpltitle = strcat(col2, 'CPT: Wo=', num2str(Wo), ', Ar=', num2str(Ar));

Z_TMet = NaN(size(X,1));
Z_Circ = NaN(size(X,1));
for k = 1:size(ModulusTable,1)
    if round(ModulusTable(k,2)) == 0 || round(ModulusTable(k,1)) == 0
        continue
    end
    i = find(round(ModulusTable(k,2)) == round(Y(:,1))); % K
    j = find(round(ModulusTable(k,1)) == round(X(1,:))); % G
    Z_TMet(i,j) = ModulusTable(k,8);
    Z_Circ(i,j) = ModulusTable(k,9);
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.1 0.2], [0.05 0.001]);
%MyMap = brewermap(21, 'Reds'); %MyMap = flipud(MyMap);
MyMap = brewermap(41, 'RdBu'); MyMap = flipud(MyMap); MyMap = MyMap(21:41,:);

% plot tensiomet data
subplot(2,5,subpl); hold on
pl_TMet = pcolor(X,Y,Z_TMet); pl_TMet.EdgeColor = 'none';
colormap(MyMap);
caxis(plot_lims); %colorbar
ax = gca; ax.TickLength = [0.03, 0.03];
set(gca,'xscale','log','yscale','log', 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', 'layer', 'top', 'TickLabelInterpreter','latex');
grid on
ax.GridColor = [0 0 0]; ax.GridAlpha = 1; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 1; ax.MinorGridLineStyle = '--';
ttl = title(TMet_subpltitle, 'interpreter', 'latex');
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
xl = xlim; yl = ylim;
plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5);  % Left & Lower Axes
plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5);  % Right & Upper Axes
if subpl == 1
    ylabel('Dilatational modulus, $K$ [mN/m]', 'interpreter', 'latex');
    set(gca, 'XTickLabel',[]);
elseif subpl == 2 || subpl == 3 || subpl == 4
    set(gca, 'XTickLabel',[], 'YTickLabel', []);
end

ax2 = axes('Position',[0.195*subpl-0.15 0.83 0.17 0.15], 'YDir', 'reverse');
plot([0,1],[0,1], 'Color','none'); hold on
if Wo == 1 && Ar == 5
    img = imread('Img_Ar5.png');
    img = flipud(img);
    image(img,'xdata',[0.25 0.75],'ydata',[0.05 0.9]);
elseif Wo == 0.1 && Ar == 5
    img = imread('Img_Ar5.png');
    img = flipud(img);
    image(img,'xdata',[0.37 0.63],'ydata',[0.05 0.5]);
elseif Wo == 1 && Ar == 1
    img = imread('Img_Ar1.png');
    img = flipud(img);
    image(img,'xdata',[0.25 0.75],'ydata',[0.05 0.45]);
elseif Wo == 0.1 && Ar == 1
    img = imread('Img_Ar1.png');
    img = flipud(img);
    image(img,'xdata',[0.37 0.63],'ydata',[0.05 0.28]);
    text(0.25, 0.6, strcat('$\sigma_{\alpha\beta}$=', num2str(Sigma), ' mN/m,   $\Delta A/A_0$=', num2str(Strain), ' \%'),...
        'interpreter', 'latex', 'fontsize', 16);
end
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'XTick',[], 'YTick',[], 'Color', 'none'); axis off


% plot circle data
subplot(2,5,subpl+5); hold on
pl_Circ = pcolor(X,Y,Z_Circ); pl_Circ.EdgeColor = 'none';
colormap(MyMap);
caxis(plot_lims); %colorbar
ax = gca; ax.TickLength = [0.03, 0.03];
set(gca,'xscale','log','yscale','log', 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', 'layer', 'top', 'TickLabelInterpreter','latex');
grid on
ax.GridColor = [0 0 0]; ax.GridAlpha = 1; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 1; ax.MinorGridLineStyle = '--';
ttl = title(Circ_subpltitle, 'interpreter', 'latex');
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
xl = xlim; yl = ylim;
plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5);  % Left & Lower Axes
plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5);  % Right & Upper Axes

if subpl == 1
    xlabel('Shear modulus, $G$ [mN/m]', 'interpreter', 'latex');
    ylabel('Dilatational modulus, $K$ [mN/m]', 'interpreter', 'latex');
elseif subpl == 2 || subpl == 3 || subpl == 4
    xlabel('Shear modulus, $G$ [mN/m]', 'interpreter', 'latex');
    set(gca, 'YTickLabel',[]);
end

hold on
if subpl == 4
    subplot(2,5,5)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    hCB.Position = [0.82, 0.475, 0.015, 0.325];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
    hCB.Ticks = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    hCB.TickLabels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]*100;
    hCB.Title.String = '$\Delta K $ (\%)';
    hCB.Title.Interpreter = 'latex';
    subplot(2,5,10)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    hCB.Position = [0.82, 0.1, 0.015, 0.323];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
    hCB.Ticks = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    hCB.TickLabels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]*100;
    hCB.Title.String = '$\Delta K $ (\%)';
    hCB.Title.Interpreter = 'latex';
end

%Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
%StrainMeasure = DataStruc.(VarList{1,1}).Params_phys.strainmeasure;
%Strain = round(100*(1 - DataStruc.(VarList{1,1}).Params_phys.frac));
%sgtitle(strcat('Relative Error in $K_{Inv}$: $\Delta A / A$=', num2str(Strain), ...
%    '\%, $\sigma_0$=', num2str(Sigma),  ' mN/m, Strain Measure=', StrainMeasure, '\hspace{3.5cm}'), ...
%    'interpreter', 'latex', 'FontSize', 18);

end