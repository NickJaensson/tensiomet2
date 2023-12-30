function [A,b] = jacobian_rhs_elastic(params_phys,vars_sol,vars_num)
    
    D = vars_num.D;
    w = vars_num.w;
    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    sigmas = vars_sol.sigmas;
    sigmap = vars_sol.sigmap;
    lams = vars_sol.lams;
    lamp = vars_sol.lamp;
    p0 = vars_sol.p0;
    N = vars_num.N;
    C = vars_num.C;

    Kmod = params_phys.Kmod;
    Gmod = params_phys.Gmod;

    % initialize some variables 
    Z = zeros(N);            % matrix filled with zeros
    IDL = [1, zeros(1,N-1)]; % line with single one and rest zeros
    ZL = zeros(1,N);         % line completely filled with zeros
    Z1 = zeros(N,1);         % column filled with zeros
    I = eye(N);              % unit matrix

    % Eq. 1-4 Knoche, p85, eq.5.7, Eq. 5-7, eq.5.8
    % determine r from psi (incl lams)
    % A11 = C*D
    % A13 = lams*sin(psi)
    % A16 = -(C*D*r)/lams
    % b1 = lams*cos(psi) - (C*D*r)
    A11 = C*D;
    A13 = diag(lams.*sin(psi));
    A16 = -C*diag((D*r)./lams);
    b1 = lams.*cos(psi)-C*D*r;

    % determine z from psi (incl lams)
    % A22 = C*D
    % A23 = -lams*cos(psi)
    % A26 = -(C*D*z)/lams
    % b2 = lams*sin(psi) - (C*D*z)
    A22 = C*D;
    A23 = diag(-lams.*cos(psi));
    A26 = -C*diag((D*z)./lams);
    b2 = lams.*sin(psi)-C*D*z;

    % determine psi from laplace law
    % A31 = -lams*(sigmap*sin(psi))/r^2
    % A32 = lams*g*rho
    % A33 = (lams*sigmap*cos(psi))/r + (C*sigmas)*D
    % A34 = (C*D*psi)
    % A35 = lams*sin(psi)/r
    % A36 = -(C*D*psi*sigmas)/lams
    % A38 = -lams
    % b3 = lams*P - lams*(sigmap*sin(psi))/r - lams*g*rho*z - (C*D*psi*sigmas)
    A31 = -diag(lams.*sigmap.*sin(psi)./(r.^2));
    A32 = diag(lams)*params_phys.deltarho*params_phys.grav;
    A33 = diag(lams.*sigmap.*cos(psi)./r)+C*diag(sigmas)*D;
    A34 = C*diag(D*psi);
    A35 = diag(lams.*sin(psi)./r);
    A36 = -C*diag((D*psi).*sigmas./lams);
    A38 = -lams;
    b3 = lams*p0 - lams.*sigmap.*sin(psi)./r ...
                  - lams.*z*params_phys.deltarho*params_phys.grav -C*sigmas.*(D*psi);

    % A81 = (2*int*lams*r*pi*sin(psi))
    % A83 = (int*lams*r^2*pi*cos(psi))
    % A86 = (int*r^2*pi*sin(psi))
    % b8 = V*C - (int*lams*r^2*pi*sin(psi))
    if params_phys.compresstype == 1
        % determine pressure - use volume
        wdef = w.*lams'/C; 
        A81 = 2*pi*wdef.*(r.*sin(psi))';
        A83 =   pi*wdef.*(r.^2.*cos(psi))';
        A86 =   pi*wdef.*((r.^2).*sin(psi))';
        b8 =   -pi*wdef*((r.^2).*sin(psi))+params_phys.volume;
    else
        % determine pressure - use area
        error('area compression not implemented')    
    end

    % boundary condition r(0) = 0
    A11(1,:) = IDL;
    A13(1,:) = ZL;
    A16(1,:) = ZL;
    b1(1) = -r(1);

    % boundary condition z(s0) = 0
    A22(1,:) = fliplr(IDL);
    A23(1,:) = ZL;
    A26(1,:) = ZL;
    b2(1) = -z(end);

    % boundary condition phi(0) = 0
    A31(1,:) = ZL;
    A32(1,:) = ZL;
    A33(1,:) = IDL;
    A34(1,:) = ZL;
    A35(1,:) = ZL;
    A36(1,:) = ZL;
    A38(1,:) = 0;
    b3(1) = -psi(1);

    % determine sigmas from projection of force balance
    % A41 = (C*D*sigmas)
    % A43 = -lams*sin(psi)*(sigmas - sigmap)
    % A44 = lams*cos(psi) + (C*r)*D
    % A45 = -lams*cos(psi)
    % A46 = -(C*r*D*sigmas)/lams
    % b4 = - lams*cos(psi)*(sigmas - sigmap) - (C*r*D*sigmas)    
    A41 = C*diag(D*sigmas);
    A43 = diag(lams.*sin(psi).*(sigmap-sigmas));
    A44 = diag(lams.*cos(psi))+C*diag(r)*D;
    A45 = -diag(lams.*cos(psi));
    A46 = -C*diag(r.*(D*sigmas)./lams);
    b4 = -C*r.*(D*sigmas)+lams.*cos(psi).*(sigmap-sigmas);

    switch params_phys.strainmeasure
    
    case 'pepicelli'

        % define some convenient variables
        lamsm1 = lams.^(-1); lamsm2 = lams.^(-2); lamsm3 = lams.^(-3);
        lampm1 = lamp.^(-1); lampm2 = lamp.^(-2); lampm3 = lamp.^(-3);
        J = lams.*lamp;
        K = params_phys.Kmod;
        G = params_phys.Gmod;
        
        % determine sigma^r
        % A54 = 1
        % A56 = (K*log(lams*lamp))/(lams^2*lamp) - K/(lams^2*lamp) - G/lams^3
        % A57 = G/lamp^3 - K/(lams*lamp^2) + (K*log(lams*lamp))/(lams*lamp^2)
        % b5 = gamma - sigmas - (G*(1/lams^2 - 1/lamp^2))/2 + (K*log(lams*lamp))/(lams*lamp)
        A54 = eye(N);
        A56 = diag(K*log(J).*lamsm2.*lampm1 - K*lamsm2.*lampm1 - G*lamsm3);
        A57 = diag(G*lampm3 - K*lamsm1.*lampm2 + K*log(J).*lamsm1.*lampm2);
        b5 = params_phys.sigma - sigmas - G*(lamsm2-lampm2)/2 + K*log(J)./J;
        
        % determine lambda^s
        % A65 = 1
        % A66 = (K*log(lams*lamp))/(lams^2*lamp) - K/(lams^2*lamp) + G/lams^3
        % A67 = - G/lamp^3 - K/(lams*lamp^2) + (K*log(lams*lamp))/(lams*lamp^2)
        % b6 = gamma - sigmap + (G*(1/lams^2 - 1/lamp^2))/2 + (K*log(lams*lamp))/(lams*lamp)
        A65 = eye(N);
        A66 = diag(K*log(J).*lamsm2.*lampm1 - K*lamsm2.*lampm1 + G*lamsm3);
        A67 = diag(-G*lampm3 - K*lamsm1.*lampm2 + K*log(J).*lamsm1.*lampm2);
        b6 = params_phys.sigma - sigmap + G*(lamsm2-lampm2)/2 + K*log(J)./J;

    end

    % A71 = 1
    % A77 = -rstar
    % b7 = -r + lamp*rstar
    % determine lambda^r
    A71 = eye(N);
    A77 = -diag(vars_sol.r_star);
    b7 = -r+lamp.*vars_sol.r_star;

    % boundary condition dsigmas/ds(0) = 0
    % NOTE: this BC is included in the Newton-Raphson iteration
    A41(1,:) = ZL;
    A43(1,:) = ZL;
    A44(1,:) = (1/lams(1))*vars_num.C*D(1,:);
    A45(1,:) = ZL;
    A46(1,:) = ZL;
    A46(1,1) = -(1/lams(1)^2)*vars_num.C*(D(1,:)*sigmas);
    b4(1) = -(1/lams(1))*vars_num.C*(D(1,:)*sigmas);

    % boundary condition lamp(s0) = 1
    A71(end,:) = ZL;
    A77(end,:) = fliplr(IDL);
    b7(end) =  1.0 - lamp(end);

    % boundary condition sigmas(0)=sigmap(0)
    A54(1,:) = IDL;
    A56(1,:) = 0;
    A57(1,:) = 0;
    A55 = Z;
    A55(1,:) = -IDL;
    b5(1) = 0;

    % boundary condition lams(0)=lamp(0)
    A71(1,:) = 0;
    A77(1,:) = IDL;
    A76 = Z;
    A76(1,:) = -IDL;
    b7(1) = 0;

    % combine matrices
    A = [[A11,   Z, A13,   Z,    Z, A16,   Z,  Z1]; ...
         [  Z, A22, A23,   Z,    Z, A26,   Z,  Z1]; ...
         [A31, A32, A33, A34,  A35, A36,   Z, A38]; ...
         [A41,   Z, A43, A44,  A45, A46,   Z,  Z1]; ...
         [  Z,   Z,   Z, A54,  A55, A56, A57,  Z1]; ...
         [  Z,   Z,   Z,   Z,  A65, A66, A67,  Z1]; ...
         [A71,   Z,   Z,   Z,    Z,   Z, A77,  Z1]; ...
         [A81, Z1', A83,  Z1', Z1', A86, Z1',   0]];

    b = [b1;b2;b3;b4;b5;b6;b7;b8];

end