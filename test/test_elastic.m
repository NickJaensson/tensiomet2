close all; clear

% load the parameter values

parameters_numerical;
parameters_simple;
parameters_elastic;

% run test for both Hencky model and Pepicelli model

for ii = 1:2

    if ii == 1
        params_phys.strainmeasure = 'pepicelli';
    elseif ii == 2    
        params_phys.strainmeasure = 'hencky';
    end

    [vars_num_ref, vars_sol_ref] = gen_single_drop(params_phys, params_num);
    
    [vars_num,vars_sol] = gen_single_drop_elastic(params_phys, params_num,...
        vars_num_ref, vars_sol_ref);
    
    [volume,area] = calculate_volume_area(vars_sol, vars_num, true);
    
    % compare to old values (gen-pendant-drop before refactoring:
    eps_test = 1e-10; 
    if strcmp(params_phys.strainmeasure, 'pepicelli')
        assert ( abs(volume-12.8000000000001) < eps_test );
        assert ( abs(area-22.5156483902096) < eps_test );
        assert ( abs(vars_sol.p0-2.02056164104124) < eps_test );
        assert ( abs(max(vars_sol.sigmas)-3.33958761227925) < eps_test );
        assert ( abs(max(vars_sol.sigmap)-3.86864491619756) < eps_test );
    elseif strcmp(params_phys.strainmeasure, 'hencky')
        assert ( abs(volume-12.8000000000161) < eps_test );
        assert ( abs(area-22.4710131205708) < eps_test );
        assert ( abs(vars_sol.p0-2.14853807101509) < eps_test );
        assert ( abs(max(vars_sol.sigmas)-3.42627115493715) < eps_test );
        assert ( abs(max(vars_sol.sigmap)-3.88525423098743) < eps_test );
    else
        error('Error test_elastic: incorrect strain measure')
    end

end

disp('All tests passed!')