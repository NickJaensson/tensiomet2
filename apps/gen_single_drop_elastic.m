% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

addpath('../src/')

gen_single_drop

% physical pa rameters
params.Kmod = 1;         % elastic dilational modulus [mN/m]
params.Gmod = 1;          % elastic shear modulus [mN/m]
params.compresstype = 1;  % 1: compress the volume    other: compress the area
params.frac = [0.9];      % compute elastic stresses for these compressions
params.strainmeasure = 'pepicelli'; % which elastic constitutive model

r0 = r; z0 = z;

params.C = C;
N=params.N;
D=params.D;
w=params.w;

% intial aread
params.area0 = pi*2*params.w*(r)/C;

% initialize the surface strains ans tresses
lamp = ones(N,1); lams = lamp;
taus = params.sigma*ones(N,1); taup = taus;

clear itervars

itervars.r0 = r0;
itervars.z0 = z0;

for ii = 1:length(params.frac)

    % determine the current target volume/area
    if params.compresstype == 1
        params.volume = params.volume0*params.frac(ii);
    else
        params.area = params.area0*params.frac(ii);
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*N+2,1);
    itervars.r = r; itervars.z = z; itervars.psi = psi;
    itervars.taus = taus; itervars.taup = taup;
    itervars.lams = lams; itervars.lamp = lamp;    
    itervars.p0 = p0; 

    while rms(u) > 1e-10

        iter = iter + 1;
    
        if iter > 1200 
            error('iter > 1200!');
        end
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_elastic(params,itervars);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        itervars.r   = itervars.r + params.alpha*u(1:N);
        itervars.z   = itervars.z + params.alpha*u(N+1:2*N);
        itervars.psi = itervars.psi + params.alpha*u(2*N+1:3*N);    
        itervars.taus = itervars.taus + params.alpha*u(3*N+1:4*N);
        itervars.taup = itervars.taup + params.alpha*u(4*N+1:5*N);
        itervars.lams = itervars.lams + params.alpha*u(5*N+1:6*N);
        itervars.lamp = itervars.lamp + params.alpha*u(6*N+1:7*N);
        itervars.p0  = itervars.p0 + params.alpha*u(end);

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

    end

    % extract the solution variables
    r = itervars.r; z = itervars.z; psi = itervars.psi;
    taus = itervars.taus; taup = itervars.taup; 
    lams = itervars.lams; lamp = itervars.lamp; 
    p0 = itervars.p0;

    % calculate the volume and the area
    volume = pi*params.w*(r.^2.*sin(psi))/C;
    area = pi*2*params.w*(r)/C;
    
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
    rmax = max([r0',rr']);
    zmin = min([z0',zz']);
    
    % rescale the plot
    xlim([0 1.2*rmax]);
    ylim([1.2*zmin 0]);
    set(gca,'DataAspectRatio',[1 1 1])
    
    % compute the curvatures (NOTE: d/ds operator is given by C*D, see Nagel)
    kappas = (C*D*psi)./lams;
    kappap = sin(psi)./r;
    kappap(1) = kappas(1);
    
    % NOTE: there are three coordinates involved: s0 (guessed domain length), 
    % s* (domain for isotropic solution, which is also the reference domain 
    % for the elastic problem) and s (domain in deformed state). The
    % grid (and thus the differentation/integration operators) is defined for 
    % the guessed domain. To obtain the other domains, we
    % use: s* = s0 / C  and  s = \int lambdas ds* = \int lambdas ds / C
    
    % construct the integration matrix from the integration vector
    wmat = repmat(w,N,1);
    wmat = tril(wmat);
    
    % compute the value of s in the deformed state
    sdef = wmat*lams/C;

end
