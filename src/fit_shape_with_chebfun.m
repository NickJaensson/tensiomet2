function [vars_sol_fit, vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise, zz_noise, params_num, cheb_num)
    % FIT_SHAPE_WITH_CHEBFUN Fits noisy shape using Chebyshev polynomials.
    %
    % This function interpolates and smooths a noisy shape using Chebyshev
    % polynomials to obtain a fitted representation.
    %
    % INPUTS:
    %   rr_noise   - Noisy radial coordinate values
    %   zz_noise   - Noisy axial coordinate values
    %   params_num - Structure containing numerical parameters
    %   cheb_num   - Optional: force number of Cheby functions
    %
    % OUTPUTS:
    %   vars_sol_fit - Structure with fitted solution variables
    %   vars_num_fit - Structure with numerical grid parameters

    if nargin < 4
        cheb_num = 0;
    end

    % get continuous s around full shape
    ds = sqrt((zz_noise(2:end)-zz_noise(1:end-1)).^2 + ...
              (rr_noise(2:end)-rr_noise(1:end-1)).^2);
    ss_noise = zeros(length(rr_noise),1);
    for i=1:length(ds)
        ss_noise(1+i) = sum(ds(1:i));
    end
    
    % find smallest stepsize and interpolate on smaller grid
    dsmin = min(diff(ss_noise));
    ssb = linspace(0,ss_noise(end),ceil(0.5*ss_noise(end)/dsmin));
    rrb = interp1(ss_noise,rr_noise,ssb,'spline');
    zzb = interp1(ss_noise,zz_noise,ssb,'spline');
    
    % fit the points using Chebyshev polynomials
    if cheb_num == 0
        frr = chebfun(rrb',[ssb(1),ssb(end)],'equi', ...
            'eps',params_num.eps_cheb);
        fzz = chebfun(zzb',[ssb(1),ssb(end)],'equi', ...
            'eps', params_num.eps_cheb);
    else
        frr = chebfun(rrb',[ssb(1),ssb(end)],cheb_num,'equi');
        fzz = chebfun(zzb',[ssb(1),ssb(end)],cheb_num,'equi');
    end

    coeffs_r = chebcoeffs(frr);
    for i = 1:size(coeffs_r,1)
        if rem(i,2) ~= 0 % odd number
            coeffs_r(i) = 0;
        end
    end
    
    coeffs_z = chebcoeffs(fzz);
    for i = 1:size(coeffs_z,1)
        if rem(i,2) == 0 % even number
            coeffs_z(i) = 0;
        end
    end
 
    if cheb_num == 0
        fr = chebfun(coeffs_r,[ssb(1),ssb(end)],'coeffs');
        fz = chebfun(coeffs_z,[ssb(1),ssb(end)],'coeffs');
    else
        fr = chebfun(coeffs_r,[ssb(1),ssb(end)],cheb_num,'coeffs');
        fz = chebfun(coeffs_z,[ssb(1),ssb(end)],cheb_num,'coeffs');
    end
  
    vars_sol_fit.r = gridsample(fr,params_num.N,[ssb(end)/2,ssb(end)]);
    vars_sol_fit.z = gridsample(fz,params_num.N,[ssb(end)/2,ssb(end)]);
    
    % we use the length of the FITTED shape for the new numerical domain
    new_length = integral(sqrt(diff(fr)^2+diff(fz)^2))/2;
    
    % now the mesh is for half of the shape (similar to forward problem)
    vars_num_fit = numerical_grid(params_num,[0,new_length]);
    dummy.C = 1;
    vars_num_fit = update_numerical_grid(dummy, vars_num_fit, false);
    
    vars_sol_fit.psi = atan2(vars_num_fit.Ds*vars_sol_fit.z, ...
                             vars_num_fit.Ds*vars_sol_fit.r);

end
