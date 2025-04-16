% dimensionfull input parameters
Kmod_guess = 2;
Gmod_guess = 20;
sigma_guess = 60;
p0_guess = 2*sigma_guess/rneedle;  

% numerical parameters for inverse problem
params_num.eps_cheb = 1e-5;      % error for describing the shape

% parameters for Young-Laplace inverse problem
params_num.eps_inv_yl = 1e-5;    % convergence critertion inverse problem
params_num.sigma_guess = sigma_guess/(deltarho*abs(grav)*rneedle^2);     % guess for interfacial tension value
params_num.p0_guess = p0_guess/(deltarho*abs(grav)*rneedle);         % guess for pressure
params_num.alpha_yl = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv_yl = 500; % maximum number of iteration steps invers

% parameters for SFE inverse problem
params_num.eps_inv_sfe = 1e-6;   % convergence critertion inverse problem
params_num.K_guess = Kmod_guess/(deltarho*abs(grav)*rneedle^2);          % guess for dilational modulus
params_num.G_guess = Gmod_guess/(deltarho*abs(grav)*rneedle^2);          % guess for shear moduls
params_num.alpha_sfe = 0.5;      % relaxation parameter in inverse problem
params_num.maxiter_inv_sfe = 5000;% maximum number of iteration steps invers