% numerical parameters for inverse problem
params_num.eps_cheb = 1e-4;      % error for describing the shape

% parameters for Young-Laplace inverse problem
params_num.eps_inv_yl = 1e-5;    % convergence critertion inverse problem
params_num.sigma_guess = 10;     % guess for interfacial tension value
params_num.p0_guess = 5;         % guess for pressure
params_num.alpha_yl = 1.0;       % relaxation parameter in inverse problem
params_num.maxiter_inv_yl = 500; % maximum number of iteration steps invers

% parameters for SFE inverse problem
params_num.eps_inv_sfe = 1e-4;   % convergence critertion inverse problem
params_num.K_guess = 1;          % guess for dilational modulus
params_num.G_guess = 1;          % guess for shear moduls
params_num.alpha_sfe = 0.5;      % relaxation parameter in inverse problem
params_num.maxiter_inv_sfe = 500;% maximum number of iteration steps invers