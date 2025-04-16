% Make Operating Windows plot


%subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.1 0.22], [0.05 0.001]);
%subplot(2,5,subpl); hold on

figure
%loglog(1,1,'w');


%patch(X,Y,C);
x_lims = [1 1 100 100]; y_lims = [1 100 100 1];
x_CPT = [1 1 100 100]; y_CPT = [100 2 20 100];
x_SFE = [2 100 100 2]; y_SFE = [100 100 10 10];
x_SFE_Avg = [2 100 100 2]; y_SFE_Avg = [100 100 2 2];

MyMap = brewermap(10, 'RdPu');
col_CPT = [37, 219, 186]/255;
col_CPT2 = [35, 153, 136]/255;
col_SFE = [219, 37, 143]/255;
col_SFE2 = [196, 31, 127]/255;
% col_CPT = 'r';
% col_CPT2 = [186, 20, 20]/255;
% col_SFE = [37, 219, 192]/255;
% col_SFE2 = [35, 153, 136]/255;
col_SFE_Avg = MyMap(10,:);
col_SFE_Avg2 = [110, 29, 133]/255;
col_edge = [200 200 200]/255;

col_CPT = 'r';
col_CPT2 = [219, 31, 31]/255;
col_SFE = 'b';
col_SFE2 = [13, 47, 148]/255;
col_SFE_Avg = [37, 219, 186]/255;
col_SFE_Avg2 = [35, 153, 136]/255;
col_edge = [200 200 200]/255;

%x = logspace(0.3,2); y = x;
%[X,Y] = meshgrid(x,y);
%Z = X;
%h = surf(X,Y,Z, 'FaceAlpha', 0.5); hold on
%set(h, 'FaceColor','none','EdgeColor','k', 'FaceAlpha', 0.5);

alpha_CPT = 0.15; 
alpha_SFE = 0.2;
alpha_SFE_Avg = 0.15;
LW = 2.5;
fill(x_CPT, y_CPT, col_CPT, 'EdgeColor', col_CPT, 'FaceColor', col_CPT, 'FaceAlpha', alpha_CPT, 'EdgeAlpha', 0.5, 'LineWidth', LW); hold on
fill(x_SFE, y_SFE, col_SFE, 'EdgeColor', col_SFE, 'FaceColor', col_SFE, 'FaceAlpha', alpha_SFE, 'EdgeAlpha', 0.5, 'LineWidth', LW);
fill(x_SFE_Avg, y_SFE_Avg, col_SFE_Avg2, 'EdgeColor', col_SFE_Avg, 'FaceColor', col_SFE_Avg, 'FaceAlpha', alpha_SFE_Avg, 'EdgeAlpha', 0.5, 'LineWidth', LW);
fill(x_lims, y_lims, 'w', 'EdgeColor', 'k', 'FaceAlpha', 0, 'LineWidth', 1); 
text(1.1, 3.5, 'CPT', 'Color', col_CPT2, 'FontSize', 16, 'fontweight', 'bold');
text(55, 12, 'SFE', 'Color', col_SFE2, 'FontSize', 16, 'fontweight', 'bold');
text(10.5, 2.5, 'SFE-Avg', 'Color', col_SFE_Avg2, 'FontSize', 16, 'fontweight', 'bold');

xl = xlim; yl = ylim;
plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5);  % Left & Lower Axes
plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5);  % Right & Upper Axes


ax = gca; ax.TickLength = [0.03, 0.03];
set(gca,'xscale','log','yscale','log', 'Xlim', [1 100], 'Ylim', [1 100], 'LineWidth', 0.5, 'FontSize',14, 'TickDir','out', 'layer', 'top', 'TickLabelInterpreter','latex');
grid on
ax.GridColor = [0 0 0]; ax.GridAlpha = 1; ax.GridLineStyle = '--';
ax.MinorGridColor = [0 0 0]; ax.MinorGridAlpha = 1; ax.MinorGridLineStyle = '--';
xlabel('Shear modulus, $G$ [mN/m]', 'interpreter', 'latex');
ylabel('Dilatational modulus, $K$ [mN/m]', 'interpreter', 'latex');

