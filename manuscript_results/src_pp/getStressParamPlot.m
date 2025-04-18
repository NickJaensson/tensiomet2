function [] = getStressParamPlot(DataStruc, VarList, StressParam_table, X, Y, subpl)

%plot_lims = [0 0.5];
plot_lims = [0 0.2];

% make exceptions for:
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.9, stress = Balemans: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Bal90')
    idx = find(StressParam_table(:,1) == 19);
    StressParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.8, stress = Balemans: 
% change G = 49 to G = 50, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal80')
    idx = find(StressParam_table(:,1) == 49);
    StressParam_table(idx,1) = 50;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Pepicelli: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Pep90')
    idx = find(StressParam_table(:,1) == 9);
    StressParam_table(idx,1) = 10;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.9, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Pep90')
    idx = find(StressParam_table(:,1) == 19);
    StressParam_table(idx,1) = 20;
end
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.8, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Pep80')
    idx = find(StressParam_table(:,1) == 19);
    StressParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.95, stress = Balemans: 
% change G = 18 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal95')
    idx = find(StressParam_table(:,1) == 18);
    StressParam_table(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Balemans: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Bal90')
    idx = find(StressParam_table(:,1) == 9);
    StressParam_table(idx,1) = 10;
end


Z = NaN(size(X,1));
for k = 1:size(StressParam_table,1)
    if round(StressParam_table(k,2)) == 0 || round(StressParam_table(k,1)) == 0
        continue
    end
    i = find(round(StressParam_table(k,2)) == round(Y(:,1))); % K
    j = find(round(StressParam_table(k,1)) == round(X(1,:))); % G
    Z(i,j) = StressParam_table(k,3);
end
%subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.2 0.2], [0.05 0.1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.1 0.22], [0.05 0.001]);
%subplot(1,5,subpl); hold on
subplot(2,5,subpl); hold on
%MyMap = brewermap(9, 'Reds'); %MyMap = flipud(MyMap);
MyMap = brewermap(21, 'PiYG'); MyMap = flipud(MyMap); MyMap = MyMap(11:21,:);
pl = pcolor(X,Y,Z); pl.EdgeColor = 'none';
%pl.FaceColor = 'interp';
colormap(MyMap);
caxis(plot_lims); %colorbar
ax = gca; ax.TickLength = [0.03, 0.03];
set(gca,'xscale','log','yscale','log', 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', 'layer', 'top', 'TickLabelInterpreter','latex');
grid on
ax.GridColor = [0 0 0]; ax.GridAlpha = 1; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 1; ax.MinorGridLineStyle = '--';
xlabel('Shear modulus, $G$ [mN/m]', 'interpreter', 'latex');
if subpl == 1
    ylabel('Dilatational modulus, $K$ [mN/m]', 'interpreter', 'latex');
else
    set(gca, 'YTickLabel',[]);
end

xl = xlim; yl = ylim;
plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5);  % Left & Lower Axes
plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5);  % Right & Upper Axes
sigma0 = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
Wo = DataStruc.(VarList{1,1}).Params_phys.Wo_paper;
Ar = DataStruc.(VarList{1,1}).Params_phys.Ar_paper;
switch subpl
    case 1
        if sigma0 == 20
            col = '(a) \,';
        elseif sigma0 == 60
            col = '(e) \,';
        else
            col = '(a) \,';
        end
    case 2
        if sigma0 == 20
            col = '(b) \,';
        elseif sigma0 == 60
            col = '(f) \,';
        else
            col = '(b) \,';
        end
    case 3
        if sigma0 == 20
            col = '(c) \,';
        elseif sigma0 == 60
            col = '(g) \,';
        else
            col = '(c) \,';
        end
    case 4
        if sigma0 == 20
            col = '(d) \,';
        elseif sigma0 == 60
            col = '(h) \,';
        else
            col = '(d) \,';
        end
end
%TitleLabel = matlab.lang.makeValidName(strcat('DAP: ', DilParam));
strainlabel = num2str(round((1 - DataStruc.(VarList{1,1}).Params_phys.frac)*100));
ttl = title(strcat({strcat(col, '\textbf{SAP:} Wo=', num2str(Wo), ', Ar=', num2str(Ar), ',') ; ...
    strcat('\qquad $\sigma_{\alpha\beta}$=', num2str(sigma0), ' mN/m, $\Delta A/A_0$=', strainlabel, ' \%')}), 'interpreter', 'latex');
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
hold on

if sigma0 == 20
ax2 = axes('Position',[0.195*subpl-0.15 0.83 0.17 0.15], 'YDir', 'reverse');
plot([0,1],[0,1], 'Color','none'); hold on
if Wo == 1 && Ar == 5
    img = imread('Img_Ar5.png');
    img = flipud(img);
    image(img,'xdata',[0.25 0.75],'ydata',[0.05 0.9]+0.12);
elseif Wo == 0.1 && Ar == 5
    img = imread('Img_Ar5.png');
    img = flipud(img);
    image(img,'xdata',[0.37 0.63],'ydata',[0.05 0.5]+0.12);
elseif Wo == 1 && Ar == 1
    img = imread('Img_Ar1.png');
    img = flipud(img);
    image(img,'xdata',[0.25 0.75],'ydata',[0.05 0.45]+0.12);
elseif Wo == 0.1 && Ar == 1
    img = imread('Img_Ar1.png');
    img = flipud(img);
    image(img,'xdata',[0.37 0.63],'ydata',[0.05 0.28]+0.12);
end
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'XTick',[], 'YTick',[], 'Color', 'none'); axis off
end

if subpl == 4
    subplot(1,4,5)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    %hCB.Position = [0.76, 0.2, 0.015, 0.61];
    hCB.Position = [0.84, 0.461, 0.015, 0.318];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
end

%Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
%StrainMeasure = DataStruc.(VarList{1,1}).Params_phys.strainmeasure;
%Strain = round(100*(1 - DataStruc.(VarList{1,1}).Params_phys.frac));
%sgtitle(strcat('Stress Anisotropy Parameter (SAP): $\Delta A / A$=', num2str(Strain), ...
%    '\%, $\sigma_0$=', num2str(Sigma),  ' mN/m, Strain Measure=', StrainMeasure, '\hspace{3.5cm}'), ...
%    'interpreter', 'latex', 'FontSize', 18);

end