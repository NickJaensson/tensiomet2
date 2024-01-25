close all; clear

syms sigma rho g V A C K G

syms r z psi sigmas sigmap lams lamp P psi_prime r_prime z_prime sigmas_prime sigmap_prime lams_prime lamp_prime P_prime
syms dr dz dpsi dsigmas dsigmap dlams dlamp dP dpsi_prime dr_prime dz_prime dsigmas_prime dsigmap_prime dlams_prime dlamp_prime dP_prime
syms kappat kappas int J rstar

kappat = sin(psi)/r;
kappas = C*psi_prime/lams;
J = lams*lamp;

% Pepicelli
% f{1} = sigmas - sigma - (K/J)*log(J) - (G/2)*(1/lamp^2 - 1/lams^2);
% f{2} = sigmap - sigma - (K/J)*log(J) + (G/2)*(1/lamp^2 - 1/lams^2);

% Hencky
f{1} = sigmas - sigma - K*log(J) - G*log(lams/lamp);
f{2} = sigmap - sigma - K*log(J) - G*log(lamp/lams);

% Hookean (linear)
% NOTE: can be derived from the Hencky model by linearizing the Hencky
% model around lams=1 and lamp=1: 
%     log(lams*lamp) ~= (lams-1)+(lamp-1)
%     log(lams/lamp) ~= (lams-1)-(lamp-1)
% f{1} = sigmas - sigma - (K+G)*(lams-1) - (K-G)*(lamp-1);
% f{2} = sigmap - sigma - (K+G)*(lamp-1) - (K-G)*(lams-1);

for i=1:length(f)

    ff = f{i};

    str = 'A'+string(i)+'1 = '+string(diff(ff,lams))  +' + '+string(diff(ff,lams_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'2 = '+string(diff(ff,K))  +' + '+string(diff(ff,z_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'3 = '+string(diff(ff,G))+' + '+string(diff(ff,psi_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'b'+string(i)+' = '+string(-ff);
    fprintf(replace_prime(str)+'\n')

    disp(' ');

end

fprintf('Note: we assume the solution vector to be:\n');
fprintf('[dlams,dK,dG]\n');

function new_str = replace_prime(str)
    new_str = strrep(str,'r_prime','D*r');
    new_str = strrep(new_str,'z_prime','D*z');
    new_str = strrep(new_str,'psi_prime','D*psi');
    new_str = strrep(new_str,'sigmas_prime','D*sigmas');
    new_str = strrep(new_str,'sigmap_prime','D*sigmap');
    new_str = strrep(new_str,'lams_prime','D*lams');
    new_str = strrep(new_str,'lamp_prime','D*lamp');
end