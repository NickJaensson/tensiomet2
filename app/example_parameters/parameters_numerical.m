% numerical parameters
params_num.N = 40;          % grid points for calculation
params_num.eps_fw = 1e-12;  % convergence criterion forward: rms(u) < eps
                            % NOTE: eps_fw used for simple and elastic
params_num.maxiter = 100;   % maximum number of iteration steps
                            % NOTE: maxiter used for simple and elastic