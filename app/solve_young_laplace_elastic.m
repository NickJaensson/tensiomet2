function [itervars,params] = solve_young_laplace_elastic(itervars,params)

    % determine the current target volume/area
    if params.compresstype == 1
        params.volume = params.volume0*params.frac;
    else
        params.area = params.area0*params.frac;
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
        [A,b] = jacobian_rhs_elastic(params,itervars);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        itervars.r   = itervars.r + params.alpha*u(1:params.N);
        itervars.z   = itervars.z + params.alpha*u(params.N+1:2*params.N);
        itervars.psi = itervars.psi + params.alpha*u(2*params.N+1:3*params.N);    
        itervars.sigmas = itervars.sigmas + params.alpha*u(3*params.N+1:4*params.N);
        itervars.sigmap = itervars.sigmap + params.alpha*u(4*params.N+1:5*params.N);
        itervars.lams = itervars.lams + params.alpha*u(5*params.N+1:6*params.N);
        itervars.lamp = itervars.lamp + params.alpha*u(6*params.N+1:7*params.N);
        itervars.p0  = itervars.p0 + params.alpha*u(end);

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

    end

    % NOTE: there are three coordinates involved: s0 (guessed domain length), 
    % s* (domain for isotropic solution, which is also the reference domain 
    % for the elastic problem) and s (domain in deformed state). The
    % grid (and thus the differentation/integration operators) is defined for 
    % the guessed domain. To obtain the other domains, we
    % use: s* = s0 / C  and  s = \int lambdas ds* = \int lambdas ds / C

    % construct the integration matrix from the integration vector
    params.wmat = repmat(params.w,params.N,1);
    params.wmat = tril(params.wmat);

    % compute the value of s in the deformed state
    params.sdef = params.wmat*itervars.lams/params.C;
    
end