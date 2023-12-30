function [r, z, s] = guess_shape(params_phys, Npoints)
% GUESS_SHAPE Predicts the initial shape of a droplet based on physical 
% parameters.
%
% This function calculates the initial guess of the shape of a droplet
% by using two different approaches depending on the Worthinton number (Wo).
% If Wo > 0.14, an empirical approach based on Nagel's work is used.
% Else, a quarter period of a cosine function with a similar volume is used.
%
% Syntax:
%   [r, z, s] = guess_shape(params_phys, Npoints)
%
% Inputs:
%   params_phys - A structure containing the physical parameters of the 
%                 droplet: Wo (Worthinton number), sigma (surface tension), 
%                 deltarho (density difference), grav (gravity), 
%                 rneedle (needle radius), and volume0 (initial volume).
%   Npoints     - The number of points to use in the discretization of the 
%                 droplet shape.
%
% Outputs:
%   r - Radial coordinate of the droplet interface.
%   z - Axial coordinate of the droplet interface.
%   s - Arc length along the droplet interface from the apex.

    % find the initial guess of the droplet shape
    if  params_phys.Wo > 0.14
    
        % predict the droplet shape using the emperical approach from Nagel
        
        % predict the maximum length of the interface (empirical Nagel)
        sigmaprime = params_phys.sigma/ ...
            (params_phys.deltarho*params_phys.grav*params_phys.rneedle^2);

        smax = sqrt(sigmaprime)*2.0/0.8701;    
    
        s = linspace(0,smax,Npoints);
    
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
        params_phys.rneedle = params_phys.rneedle; 
        r = linspace(0,params_phys.rneedle,Npoints);
        z = -sqrt(2*params_phys.volume0/ ...
            (pi*params_phys.rneedle))*cos(pi*r/(2*params_phys.rneedle));
    
        % determine the curve length
        ds = zeros(size(r));
        ds(2:end) = sqrt((r(2:end)-r(1:Npoints-1)).^2 + ...
                               (z(2:end)-z(1:Npoints-1)).^2);

        % determine the curve coordinate
        s = cumsum(ds);
        
    end

end