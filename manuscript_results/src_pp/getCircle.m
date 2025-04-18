function [r_circlefit, z_circlefit, R_circlefit] = getCircle(r_coords, z_coords, Lfrac)

% fit circle
r = r_coords(Lfrac : length(r_coords)-Lfrac);
z = z_coords(Lfrac : length(z_coords)-Lfrac);
[rc_fit,zc_fit,R_fit,~] = circfit(r,z);



% get coordinates of circle fit
L = sqrt( (abs(r(end)-r(length(r)/2+1)))^2 + (abs(z(end)-z(length(z)/2+1)))^2 );
theta = 2*asin(L/(2*R_fit));
theta_vec = linspace(-pi/2 - theta, -pi/2 + theta, 80);
r_circlefit = R_fit*cos(theta_vec)+rc_fit;
z_circlefit = R_fit*sin(theta_vec)+zc_fit;
R_circlefit = R_fit;


end