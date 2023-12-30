function [vars_sol,params] = solve_forward_young_laplace_elastic(vars_sol,params)

    % determine the current target volume/area
    if params.compresstype == 1
        params.volume = params.volume0*params.frac;
    else
        error('area compression not implemented')    
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*params.N+2,1);

    % start the Newton-Raphson iteration
    while rms(u) > params.eps
    
        iter = iter + 1;
        
        if iter > params.maxiter
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_elastic(params,vars_sol);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        vars_sol.r   = vars_sol.r + params.alpha*u(1:params.N);
        vars_sol.z   = vars_sol.z + params.alpha*u(params.N+1:2*params.N);
        vars_sol.psi = vars_sol.psi + params.alpha*u(2*params.N+1:3*params.N);    
        vars_sol.sigmas = vars_sol.sigmas + params.alpha*u(3*params.N+1:4*params.N);
        vars_sol.sigmap = vars_sol.sigmap + params.alpha*u(4*params.N+1:5*params.N);
        vars_sol.lams = vars_sol.lams + params.alpha*u(5*params.N+1:6*params.N);
        vars_sol.lamp = vars_sol.lamp + params.alpha*u(6*params.N+1:7*params.N);
        vars_sol.p0  = vars_sol.p0 + params.alpha*u(end);

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

    end

    % the integration and differentation matrices in the deformed state
    % NOTE: this construction of Ddef is simlar to first applying D*f/C,
    % and then dividing the components by the components of (1/lams)
    params.wdef = params.w.*vars_sol.lams'/params.C; 
    params.Ddef = params.C*params.D.*repelem((1./vars_sol.lams)',params.N,1); 

    % construct the integration matrix from the integration vector
    params.wmat = repmat(params.w,params.N,1);
    params.wmat = tril(params.wmat);

    % compute the value of s in the deformed state
    params.sdef = params.wmat*vars_sol.lams/params.C;
   
end