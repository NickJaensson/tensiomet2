function [d,dd,w,s]=numerical_grid(N,domain)
    % the following functions are defined in the chebfun package. Please
    % make sure that set_paths.m points to you chebfun installation
    d = diffmat(N,1,domain);
    dd = diffmat(N,2,domain);
    w = introw(N,domain);
    s = chebpts(N,domain);
end