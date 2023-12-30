function [r,z,s] = guess_shape(params_phys,N_guess)

    % find the initial guess of the droplet shape
    if  params_phys.Wo > 0.14
    
        % predict the droplet shape using the emperical approach from Nagel
        
        % predict the maximum length of the interface (empirical Nagel)
        sigmaprime = params_phys.sigma/(params_phys.deltarho*params_phys.grav*params_phys.rneedle^2);

        smax = sqrt(sigmaprime)*2.0/0.8701;    
    
        s = linspace(0,smax,N_guess);
    
        % predict the shape of the interface (empirical Nagel)
        z = -4/3*smax/pi*(cos(pi*3/4*s/smax));
        z = z - max(z);
        r = 4/3*smax/pi*(sin(pi*3/4*s/smax));

        % scale the shape to match the radius
        scale = params_phys.rneedle/r(end);
        r = scale*r;
        z = scale*z;
    
    else
        
        % predict the droplet shape using a quarter of a period of a cosine
        % with similar volume as imposed
    
        % find the initial guess of the droplet shape
        N_guess = 1000;
        params_phys.rneedle = params_phys.rneedle; 
        r = linspace(0,params_phys.rneedle,N_guess);
        z = -sqrt(2*params_phys.volume0/(pi*params_phys.rneedle))*cos(pi*r/(2*params_phys.rneedle));
    
        % determine the curve length
        ds = zeros(size(r));
        ds(2:end) = sqrt((r(2:end)-r(1:N_guess-1)).^2 + ...
                               (z(2:end)-z(1:N_guess-1)).^2);

        % determine the curve coordinate
        s = cumsum(ds);
        
    end

end