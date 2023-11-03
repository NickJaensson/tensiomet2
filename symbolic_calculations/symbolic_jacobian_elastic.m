close all; clear

syms gamma rho g V C K G

syms r z psi sigmas sigmat lams lamt P psi_prime r_prime z_prime sigmas_prime sigmat_prime lams_prime lamt_prime P_prime
syms dr dz dpsi dsigmas dsigmat dlams dlamt dP dpsi_prime dr_prime dz_prime dsigmas_prime dsigmat_prime dlams_prime dlamt_prime dP_prime
syms kappat kappas int J rstar

kappat = sin(psi)/r;
kappas = C*psi_prime/lams;
J = lams*lamt;

f{1} = C*r_prime/lams - cos(psi);
f{2} = C*z_prime/lams - sin(psi);
f{3} = (kappat*sigmat+kappas*sigmas) - P + rho*g*z;
f{4} = r*C*sigmas_prime/lams - cos(psi)*(sigmat-sigmas);

f{5} = sigmas - gamma - (K/J)*log(J) - (G/2)*(1/lamt^2 - 1/lams^2);
f{6} = sigmat - gamma - (K/J)*log(J) + (G/2)*(1/lamt^2 - 1/lams^2);
f{7} = lamt - r/rstar;
f{8} = int*(pi*r^2*sin(psi)) - C*V/lams;

for i=1:length(f)

    ff = f{i};

    str = 'A'+string(i)+'1 = '+string(diff(ff,r))  +' + '+string(diff(ff,r_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'2 = '+string(diff(ff,z))  +' + '+string(diff(ff,z_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'3 = '+string(diff(ff,psi))+' + '+string(diff(ff,psi_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'4 = '+string(diff(ff,sigmas))+' + '+string(diff(ff,sigmas_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'5 = '+string(diff(ff,sigmat))+' + '+string(diff(ff,sigmat_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'6 = '+string(diff(ff,lams))+' + '+string(diff(ff,lams_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'7 = '+string(diff(ff,lamt))+' + '+string(diff(ff,lamt_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'8 = '+string(diff(ff,P))  +' + '+string(diff(ff,P_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'b'+string(i)+' = '+string(-ff);
    fprintf(replace_prime(str)+'\n')

    disp(' ');

end

fprintf('Note: we assume the solution vector to be:\n');
fprintf('[dr,dz,dpsi,dsigmas,dsigmat,dlams,dlamt,dP]\n');

function new_str = replace_prime(str)
    new_str = strrep(str,'r_prime','D*r');
    new_str = strrep(new_str,'z_prime','D*z');
    new_str = strrep(new_str,'psi_prime','D*psi');
    new_str = strrep(new_str,'sigmas_prime','D*sigmas');
    new_str = strrep(new_str,'sigmat_prime','D*sigmat');
    new_str = strrep(new_str,'lams_prime','D*lams');
    new_str = strrep(new_str,'lamt_prime','D*lamt');
end