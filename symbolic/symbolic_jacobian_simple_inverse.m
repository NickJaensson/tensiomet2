close all; clear

syms sigma rho g V

syms r z psi C P psi_prime r_prime  z_prime C_prime P_prime r_exp
syms dr dz dpsi dC dP dpsi_prime dr_prime dz_prime dC_prime dP_prime
syms kappat kappas int

kappat = sin(psi)/r;
kappas = psi_prime;

f{1} = C*r_prime - cos(psi);
f{2} = C*z_prime - sin(psi);
f{3} = -P+rho*g*z+sigma*(C*kappas+kappat);

for i=1:3

    ff = f{i};

    str = 'A'+string(i)+'1 = '+string(diff(ff,r))  +' + '+string(diff(ff,r_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'2 = '+string(diff(ff,z))  +' + '+string(diff(ff,z_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'3 = '+string(diff(ff,psi))+' + '+string(diff(ff,psi_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'4 = '+string(diff(ff,sigma))  +' + '+string(diff(ff,C_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    str = 'A'+string(i)+'5 = '+string(diff(ff,P))  +' + '+string(diff(ff,P_prime))+'*D';
    fprintf(replace_prime(str)+'\n')

    disp(' ');

end

fprintf('Note: we assume the solution vector to be: [dr,dz,dpsi,dsigma,dP]\n');

function new_str = replace_prime(str)
    new_str = strrep(str,'r_prime','D*r');
    new_str = strrep(new_str,'z_prime','D*z');
    new_str = strrep(new_str,'psi_prime','D*psi');
end
