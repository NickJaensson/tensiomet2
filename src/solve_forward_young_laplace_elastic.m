function [vars_sol,params_num] = solve_forward_young_laplace_elastic(vars_sol,params_phys,params_num)

    % determine the current target volume/area
    if params_phys.compresstype == 1
        params_phys.volume = params_phys.volume0*params_phys.frac;
    else
        error('area compression not implemented')    
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*params_num.N+2,1);

    % start the Newton-Raphson iteration
    while rms(u) > params_num.eps
    
        iter = iter + 1;
        
        if iter > params_num.maxiter
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_elastic(params_phys,params_num,vars_sol);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        vars_sol.r   = vars_sol.r + params_num.alpha*u(1:params_num.N);
        vars_sol.z   = vars_sol.z + params_num.alpha*u(params_num.N+1:2*params_num.N);
        vars_sol.psi = vars_sol.psi + params_num.alpha*u(2*params_num.N+1:3*params_num.N);    
        vars_sol.sigmas = vars_sol.sigmas + params_num.alpha*u(3*params_num.N+1:4*params_num.N);
        vars_sol.sigmap = vars_sol.sigmap + params_num.alpha*u(4*params_num.N+1:5*params_num.N);
        vars_sol.lams = vars_sol.lams + params_num.alpha*u(5*params_num.N+1:6*params_num.N);
        vars_sol.lamp = vars_sol.lamp + params_num.alpha*u(6*params_num.N+1:7*params_num.N);
        vars_sol.p0  = vars_sol.p0 + params_num.alpha*u(end);

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

    end

    % the integration and differentation matrices in the deformed state
    % NOTE: this construction of Ddef is simlar to first applying D*f/C,
    % and then dividing the components by the components of (1/lams)
    params_num.wdef = params_num.w.*vars_sol.lams'/params_num.C; 
    params_num.Ddef = params_num.C*params_num.D.*repelem((1./vars_sol.lams)',params_num.N,1); 

    % construct the integration matrix from the integration vector
    params_num.wmat = repmat(params_num.w,params_num.N,1);
    params_num.wmat = tril(params_num.wmat);

    % compute the value of s in the deformed state
    params_num.sdef = params_num.wmat*vars_sol.lams/params_num.C;
   
end