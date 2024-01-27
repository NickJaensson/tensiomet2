% numerical parameters for inverse problem
params_num.eps_cheb = 1e-4;   % error for describing the shape
params_num.eps_inv = 1e-5;    % convergence critertion inverse problem

% initial guesses for Young-Laplace inverse problem
params_num.sigma_guess = 10;  % guess for interfacial tension value
params_num.p0_guess = 5;      % guess for pressure

% initial guesses for SFE inverse problem
params_num.K_guess = 1;      % guess for dilational modulus
params_num.G_guess = 1;      % guess for shear moduls

params_num.alpha = 1.0;       % relaxation parameter in inverse problem
                              % NOTE: used for both YL and SFE
params_num.maxiter_inv = 500; % maximum number of iteration steps inverse
                              % NOTE: used for both YL and SFE
