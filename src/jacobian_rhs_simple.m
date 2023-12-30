function [A,b] = jacobian_rhs_simple(params_phys,params_num,vars_sol)
    
    D = params_num.D;
    w = params_num.w;
    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;
    
    % initialize some variables 
    Z = zeros(params_num.N);            % matrix filled with zeros
    IDL = [1, zeros(1,params_num.N-1)]; % line with single one and rest zeros
    ZL = zeros(1,params_num.N);         % line completely filled with zeros 
    b = ones(3*params_num.N+2,1); % solution vector and right hand side
    
    % determine r from psi
    A11 = C*D; A13 = diag(sin(psi)); A14 = D*r; b1 = -(C*D*r-cos(psi));
    
    % determine z from psi 
    A22 = C*D; A23 = diag(-cos(psi)); A24 = D*z; b2 = -(C*D*z-sin(psi));
    
    % determine psi from Laplace law
    A31 = -params_phys.sigma*diag(sin(psi)./r.^2);
    A32 = params_phys.grav*params_phys.deltarho*diag(ones(params_num.N,1));
    A33 = C*params_phys.sigma*D + params_phys.sigma*diag(cos(psi)./r);
    A34 = params_phys.sigma*(D*psi);
    A35 = -ones(params_num.N,1);
    b3 = p0-params_phys.grav*params_phys.deltarho*z-params_phys.sigma*(C*D*psi+sin(psi)./r);
    
    % impose the needle radius as a BC (imposes the domain length)
    A41 = fliplr(IDL); b4 = (params_phys.rneedle-r(end));
    
    % determine pressure - use volume
    A51 = pi*2*w.*r'.*sin(psi');
    A53 = pi*w.*r'.^2.*cos(psi');
    A54 = -params_phys.volume0;
    b5 = -(pi*w*(r.^2.*sin(psi))-C*params_phys.volume0);
    
    % boundary condition r(0) = 0
    A11(1,:) = IDL; 
    A13(1,:) = ZL; 
    A14(1) = 0;
    b1(1) = -r(1);
    
    % boundary condition z(s0) = 0
    A22(1,:) = fliplr(IDL); 
    A23(1,:) = ZL; 
    A24(1) = 0;
    b2(1) = -z(end);
    
    % boundary condition phi(0) = 0
    A31(1,:) = ZL; 
    A32(1,:) = ZL; 
    A33(1,:) = IDL; 
    A34(1,:) = 0; 
    A35(1,:) = 0;
    b3(1) = -psi(1);
    
    % assemble matrices
    Z1 = zeros(params_num.N,1);
     
    A = [[A11,   Z, A13, A14,  Z1];
       [  Z, A22, A23, A24,  Z1];
       [A31, A32, A33, A34, A35];
       [A41,  ZL,  ZL,   0,   0];
       [A51, Z1', A53, A54,   0]];
     
    b = [b1;b2;b3;b4;b5];


end