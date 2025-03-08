% plot the shapes for different Bond numbers and volumes. If the shape 
% cannot be found, nothing is displayed. 
% Blue shapes: Wo > 0.5, red shapes: Wo <= 0.5

close all; clear

% numerical parameters
params_num.N = 40;                  % grid points for calculation
params_num.eps_fw_simple = 1e-12;   % convergence criterion forward
params_num.maxiter_simple = 100;    % maximum number of iteration steps

% physical parameters for the simple droplet problem
params_phys.sigma = 1;      % surface tension
params_phys.grav = 1;       % gravitational acceleration
params_phys.rneedle = 1;    % radius of the needle

Bond_all = 1.0:-0.1:0.10;  % Bond number
Nu_all   = 2:2:22;          % dimensionless volume (Nu)

counter = 0;

for iii = 1:length(Bond_all)
    
    for jjj = 1:length(Nu_all)
                    
        counter = counter + 1;
              
        params_phys.volume0 = Nu_all(jjj);   % prescribed volume
        params_phys.deltarho = Bond_all(iii);  % density difference
        
        % calculate and display the Worthing number
        Wo = params_phys.deltarho*params_phys.grav*params_phys.volume0/...
            (2*pi*params_phys.sigma*params_phys.rneedle);
        
        disp(['Worthington number = ',num2str(Wo)]);
        
        params_phys.Wo = Wo;
        
        Bo = params_phys.grav*params_phys.rneedle^2/params_phys.sigma;
        disp(['Bond number = ',num2str(Bo)]);
        
        Nu = params_phys.volume0/params_phys.rneedle^3;
        disp(['Vol number = ',num2str(Nu)]);
        
        crash = 0;
        try
            [vars_num, vars_sol, params_phys] = gen_single_drop(params_phys, ...
                params_num, false);
        catch
            crash = 1;
        end
        
        % plot the final shape
        figure(1); subplot(length(Bond_all),length(Nu_all),counter);
        
        if ~crash
            if Wo > 0.5
              fill([vars_sol.r' 0],[-vars_sol.z' 0],'b'); 
            else
              fill([vars_sol.r' 0],[-vars_sol.z' 0],'r'); 
            end
        end
        
        xlim([0 3]);
        ylim([0 5]);
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
    
        if iii == length(Bond_all)
            text(min(xlim), min(ylim)-2.0, sprintf('%1.1f', Nu_all(jjj)),'FontSize',12);
        end

        if jjj == 1
            text(min(xlim)-4.0, min(ylim), sprintf('%1.1f', Bond_all(iii)),'FontSize',12);
        end

    end

end

ax = findobj(1,'Type','Axes');
xlabel(ax(1),'---> V/R^3','FontSize',24)
ylabel(ax(end),'---> Bo','FontSize',24)