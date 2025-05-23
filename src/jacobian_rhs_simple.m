function [A, b] = jacobian_rhs_simple(params_phys, vars_sol, vars_num)
    % JACOBIAN_RHS_SIMPLE Computes the Jacobian matrix and RHS vector for a 
    % simple shape problem.
    %
    % INPUTS:
    %   params_phys - Structure with physical parameters
    %   vars_sol    - Structure with solution variables
    %   vars_num    - Structure with numerical variables
    %
    % OUTPUTS:
    %   A - Jacobian matrix for the simple shape system
    %   b - Right-hand side vector for the simple shape system
    
    D = vars_num.D0;
    w = vars_num.w0;
    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;
    
    % initialize some variables 
    Z = zeros(vars_num.N);            % matrix filled with zeros
    IDL = [1, zeros(1,vars_num.N-1)]; % line with single one and rest zeros
    ZL = zeros(1,vars_num.N);         % line completely filled with zeros 
    b = ones(3*vars_num.N+2,1); % solution vector and right hand side
    
    % determine r from psi
    A11 = C*D; A13 = diag(sin(psi)); A14 = D*r; b1 = -(C*D*r-cos(psi));
    
    % determine z from psi 
    A22 = C*D; A23 = diag(-cos(psi)); A24 = D*z; b2 = -(C*D*z-sin(psi));
    
    % determine psi from Laplace law
    a = params_phys.a;
    b = params_phys.b;
    pmin = params_phys.pmin;
    pmax = params_phys.pmax;

    A31 = -params_phys.sigma*diag(sin(psi)./r.^2) ...
          - diag((b*exp(a - b*r)*(pmax - pmin))./(exp(a - b*r) + 1).^2);
          % +diag((a*r.*exp(-r.^2/(2*b^2)))/b^2);
    A32 = params_phys.grav*params_phys.deltarho*diag(ones(vars_num.N,1));
    A33 = C*params_phys.sigma*D + params_phys.sigma*diag(cos(psi)./r);
    A34 = params_phys.sigma*(D*psi);
    A35 = -ones(vars_num.N,1);
    b3 = p0-params_phys.grav*params_phys.deltarho*z...
        -params_phys.sigma*(C*D*psi+sin(psi)./r) ...
        + pmin + (pmax - pmin)./(exp(a - b*r) + 1);
        % +a*exp(-r.^2/(2*b^2));
    
    % impose the needle radius as a BC (imposes the domain length)
    A41 = fliplr(IDL); b4 = (params_phys.rneedle-r(end));
    
    % determine pressure - use volume
    A51 = 2*pi*w.*z'.*cos(psi');
    A52 = 2*pi*w.*r'.*cos(psi');
    A53 = -2*pi*w.*r'.*z'.*sin(psi');
    A54 = -params_phys.volume0;
    b5 = -(2*pi*w*(r.*z.*cos(psi))-C*params_phys.volume0);
    
    % boundary condition r(0) = 0
    A11(1,:) = IDL; 
    A13(1,:) = ZL; 
    A14(1) = 0;
    b1(1) = -r(1);
    
    if params_phys.impose_contact_angle    
        % boundary condition psi(s0) = contact_angle
        A31(2,:) = ZL; 
        A32(2,:) = ZL; 
        A33(2,:) = fliplr(IDL); 
        A34(2,:) = 0; 
        A35(2,:) = 0;
        b3(2) = params_phys.contact_angle-psi(end);
    else
        % boundary condition z(s0) = 0
        A22(1,:) = fliplr(IDL); 
        A23(1,:) = ZL; 
        A24(1) = 0;
        b2(1) = -z(end);
    end

    % boundary condition phi(0) = 0
    A31(1,:) = ZL; 
    A32(1,:) = ZL; 
    A33(1,:) = IDL; 
    A34(1,:) = 0; 
    A35(1,:) = 0;
    b3(1) = -psi(1);
    
    % assemble matrices
    Z1 = zeros(vars_num.N,1);
     
    A = [[A11,   Z, A13, A14,  Z1];
       [  Z, A22, A23, A24,  Z1];
       [A31, A32, A33, A34, A35];
       [A41,  ZL,  ZL,   0,   0];
       [A51, A52, A53, A54,   0]];
     
    b = [b1;b2;b3;b4;b5];


end
