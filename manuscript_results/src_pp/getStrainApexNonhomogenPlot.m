function [] = getStrainApexNonhomogenPlot(DataStruc, VarList, ApexStrainTable, X, Y, subpl)

% ApexStrainTable: | G | K | ApexStrain_fwd | ApexStrain_TMet | ApexStrain_Circ | StrainParam_TMet | StrainParam_Circ |

plot_lims = [-0.5 0.5];
strainlabel = num2str(round((1 - DataStruc.(VarList{1,1}).Params_phys.frac)*100));

switch subpl
    case 1
        col1 = '(a) \,';
        col2 = '(d) \,';
    case 2
        col1 = '(b) \,';
        col2 = '(e) \,';
    case 3
        col1 = '(c) \,';
        col2 = '(f) \,';
    otherwise
        col1 = ' ';
        col2 = ' ';
end
TMet_subpltitle = strcat(col1, 'CMD-Tensiomet: Strain=', strainlabel, '\%');
Circ_subpltitle = strcat(col2, 'CPT: Strain=', strainlabel, '\%');


% make exceptions for:
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.9, stress = Balemans: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Bal90')
    idx = find(ApexStrainTable(:,1) == 19);
    ApexStrainTable(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.8, stress = Balemans: 
% change G = 49 to G = 50, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal80')
    idx = find(ApexStrainTable(:,1) == 49);
    ApexStrainTable(idx,1) = 50;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Pepicelli: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Pep90')
    idx = find(ApexStrainTable(:,1) == 9);
    ApexStrainTable(idx,1) = 10;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.9, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Pep90')
    idx = find(ApexStrainTable(:,1) == 19);
    ApexStrainTable(idx,1) = 20;
end
% Wo = 0.1, Ar = 5, sigma = 60, strain = 0.8, stress = Pepicelli: 
% change G = 19 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo01_Ar5_Sig60_Pep80')
    idx = find(ApexStrainTable(:,1) == 19);
    ApexStrainTable(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 60, strain = 0.95, stress = Balemans: 
% change G = 18 to G = 20, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig60_Bal95')
    idx = find(ApexStrainTable(:,1) == 18);
    ApexStrainTable(idx,1) = 20;
end
% Wo = 1, Ar = 5, sigma = 20, strain = 0.9, stress = Balemans: 
% change G = 9 to G = 10, for plotting purposes
if strcmp(inputname(1), 'Wo1_Ar5_Sig20_Bal90')
    idx = find(ApexStrainTable(:,1) == 9);
    ApexStrainTable(idx,1) = 10;
end

Z_TMet = NaN(size(X,1));
Z_Circ = NaN(size(X,1));
for k = 1:size(ApexStrainTable,1)
    if round(ApexStrainTable(k,2)) == 0 || round(ApexStrainTable(k,1)) == 0
        continue
    end
    i = find(round(ApexStrainTable(k,2)) == round(Y(:,1))); % K
    j = find(round(ApexStrainTable(k,1)) == round(X(1,:))); % G
    Z_TMet(i,j) = ApexStrainTable(k,6);
    Z_Circ(i,j) = ApexStrainTable(k,7);
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.1 0.12], [0.1 0.1]);
MyMap = brewermap(11, 'BrBG'); %MyMap = flipud(MyMap);

% plot tensiomet data
subplot(2,4,subpl); hold on
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
    ylabel('Dilatational modulus, K [mN/m]', 'interpreter', 'latex');
    set(gca, 'XTickLabel',[]);
elseif subpl == 2 || subpl == 3
    set(gca, 'XTickLabel',[], 'YTickLabel', []);
end


% plot circle data
subplot(2,4,subpl+4); hold on
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
    xlabel('Shear modulus, G [mN/m]', 'interpreter', 'latex');
    ylabel('Dilatational modulus, K [mN/m]', 'interpreter', 'latex');
elseif subpl == 2 || subpl == 3
    xlabel('Shear modulus, G [mN/m]', 'interpreter', 'latex');
    set(gca, 'YTickLabel',[]);
end

hold on
if subpl == 3
    subplot(2,4,4)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    hCB.Position = [0.715, 0.525, 0.015, 0.35];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
    subplot(2,4,8)
    hCB = colorbar('West');
    caxis(plot_lims);
    set(gca,'Visible',false);
    hCB.LineWidth = 1.5;
    hCB.FontSize = 14;
    hCB.Position = [0.715, 0.11, 0.015, 0.35];
    hCB.AxisLocation = 'in';
    hCB.TickLabelInterpreter = 'latex';
end

Wo = DataStruc.(VarList{1,1}).Params_phys.Wo_paper;
Ar = DataStruc.(VarList{1,1}).Params_phys.Ar_paper;
Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
StrainMeasure = DataStruc.(VarList{1,1}).Params_phys.strainmeasure;
sgtitle(strcat('Apical Strain Nonhomogeneity: Wo=', num2str(Wo), ', Ar=', num2str(Ar), ...
    ', $\sigma_0$=', num2str(Sigma), ', StrainM=', StrainMeasure, '\hspace{3.5cm}'), ...
    'interpreter', 'latex', 'FontSize', 18);

end