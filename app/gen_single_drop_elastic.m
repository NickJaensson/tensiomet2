% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

addpath('../src/')

% generate reference shape (parameters are taken from gen_single_drop.m)
gen_single_drop

% physical parameters for the elastic problem
params.Kmod = 3;         % elastic dilational modulus
params.Gmod = 2;          % elastic shear modulus
params.compresstype = 1;  % 1: compress the volume other: compress the area
params.frac = [0.8];      % compute elastic stresses for these compressions
params.strainmeasure = 'pepicelli'; % which elastic constitutive model

params.maxiter = 2000; % OVERWRITE SINCE NR ITER DOES NOT WORK YET!
params.eps = 1e-12; % OVERWRITE SINCE NR ITER DOES NOT WORK YET!

% initialize the surface strains ans tresses
lamp = ones(params.N,1); lams = lamp;
sigmas = params.sigma*ones(params.N,1); sigmap = sigmas;

% clear itervars variable from the simple interface problem
clear itervars

% store the coordinates of the reference shape
itervars.r0 = r; itervars.z0 = z;

for ii = 1:length(params.frac)

    % determine the current target volume/area
    if params.compresstype == 1
        params.volume = params.volume0*params.frac(ii);
    else
        params.area = params.area0*params.frac(ii);
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*params.N+2,1);
    itervars.r = r; itervars.z = z; itervars.psi = psi;
    itervars.sigmas = sigmas; itervars.sigmap = sigmap;
    itervars.lams = lams; itervars.lamp = lamp;    
    itervars.p0 = p0; 

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

    % extract the solution variables
    r = itervars.r; z = itervars.z; psi = itervars.psi;
    sigmas = itervars.sigmas; sigmap = itervars.sigmap; 
    lams = itervars.lams; lamp = itervars.lamp; 
    p0 = itervars.p0;

    % calculate the volume and the area
    wdef = params.w.*lams'/C; 
    volume = pi*wdef*(r.^2.*sin(psi));
    area = pi*2*wdef*(r);
    
    disp(['volume = ', num2str(volume,15)]);
    disp(['area = ', num2str(area,15)]);
    disp(['pressure = ', num2str(p0,15)]);

    % interpolate the numerical solutions on a finer grid. 
    % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
    % though all the points (see book of Trefethen on Spectral Methods in 
    % Matlab, page  63). For plotting purposes we use a simpler interpolation 
    ss = linspace(s(1),s(end),params.Nplot)';
    rr = interp1(s,r,ss,'pchip');
    zz = interp1(s,z,ss,'pchip');
    
    % plot the droplet shape
    plot(rr,zz); 
    rmax = max([itervars.r0',rr']);
    zmin = min([itervars.z0',zz']);
    
    % rescale the plot
    xlim([0 1.2*rmax]);
    ylim([1.2*zmin 0]);
    set(gca,'DataAspectRatio',[1 1 1])
    
    % compute the curvatures (NOTE: d/ds operator is given by C*D, see Nagel)
    kappas = (C*params.D*psi)./lams;
    kappap = sin(psi)./r;
    kappap(1) = kappas(1);
    
    % NOTE: there are three coordinates involved: s0 (guessed domain length), 
    % s* (domain for isotropic solution, which is also the reference domain 
    % for the elastic problem) and s (domain in deformed state). The
    % grid (and thus the differentation/integration operators) is defined for 
    % the guessed domain. To obtain the other domains, we
    % use: s* = s0 / C  and  s = \int lambdas ds* = \int lambdas ds / C
    
    % construct the integration matrix from the integration vector
    wmat = repmat(params.w,params.N,1);
    wmat = tril(wmat);
    
    % compute the value of s in the deformed state
    sdef = wmat*lams/C;

    % figure; hold on
    % plot(sdef,itervars.sigmas)
    % plot(sdef,itervars.sigmap)
    % 
    % disp(['max(sigmas) = ', num2str(max(sigmas),15)]);
    % disp(['max(sigmap) = ', num2str(max(sigmap),15)]);

end