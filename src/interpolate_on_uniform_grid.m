function varargout = interpolate_on_uniform_grid(vars_num,f,Ngrid)

    % interpolate the numerical solutions on a uniform grid.
    % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
    % though all the points (see book of Trefethen on Spectral Methods in 
    % Matlab, page  63). For plotting purposes we use a simpler interpolation 

    s_uniform = linspace(vars_num.s(1),vars_num.s(end),Ngrid)';

    varargout{1} = s_uniform;

    for i = 1:size(f,2)
        varargout{i+1} = interp1(vars_num.s,f(:,i),s_uniform,'pchip');
    end

end