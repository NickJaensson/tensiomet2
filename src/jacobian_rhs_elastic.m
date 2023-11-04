function [A,b] = jacobian_rhs_simple(params,itervars)
    
    D = params.D;
    w = params.w;
    r = itervars.r;
    z = itervars.z;
    psi = itervars.psi;
    sigmas = itervars.taus;
    sigmat = itervars.taup;
    lams = itervars.lams;
    lamt = itervars.lamp;
    p0 = itervars.p0;
    N = params.N;
    C = params.C;

    Kmod = params.Kmod;
    Gmod = params.Gmod;

    % initialize some variables 
    Z = zeros(N);            % matrix filled with zeros
    IDL = [1, zeros(1,N-1)]; % line with single one and rest zeros
    ZL = zeros(1,N);         % line completely filled with zeros
    Z1 = zeros(N,1);
    I = eye(N);

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
    % A31 = -lams*(sigmat*sin(psi))/r^2
    % A32 = lams*g*rho
    % A33 = (lams*sigmat*cos(psi))/r + (C*sigmas)*D
    % A34 = (C*D*psi)
    % A35 = lams*sin(psi)/r
    % A36 = -(C*D*psi*sigmas)/lams
    % A38 = -lams
    % b3 = lams*P - lams*(sigmat*sin(psi))/r - lams*g*rho*z - (C*D*psi*sigmas)
    A31 = -diag(lams.*sigmat.*sin(psi)./(r.^2));
    A32 = diag(lams)*params.deltarho*params.grav;
    A33 = diag(lams.*sigmat.*cos(psi)./r)+C*diag(sigmas)*D;
    A34 = C*diag(D*psi);
    A35 = diag(lams.*sin(psi)./r);
    A36 = -C*diag((D*psi).*sigmas./lams);
    A38 = -lams;
    b3 = lams*p0 - lams.*sigmat.*sin(psi)./r ...
                  - lams.*z*params.deltarho*params.grav -C*sigmas.*(D*psi);

    % A81 = (2*int*lams*r*pi*sin(psi))
    % A83 = (int*lams*r^2*pi*cos(psi))
    % A86 = (int*r^2*pi*sin(psi))
    % b8 = V*C - (int*lams*r^2*pi*sin(psi))
    if params.compresstype == 1 || ii == 1
        % determine pressure - use volume
        wdef = w.*lams'/C; 
        A81 = 2*pi*wdef.*(r.*sin(psi))';
        A83 =   pi*wdef.*(r.^2.*cos(psi))';
        A86 =   pi*wdef.*((r.^2).*sin(psi))';
        b8 =   -pi*wdef*((r.^2).*sin(psi))+params.volume;
    else
        % determine pressure - use area
        A81 = 2*w.*lams';
        A83 = zeros(1,N);
        b8 = -(2*w*(r.*lams)-C*params.area/pi);
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
    % A43 = -lams*sin(psi)*(sigmas - sigmat)
    % A44 = lams*cos(psi) + (C*r)*D
    % A45 = -lams*cos(psi)
    % A46 = -(C*r*D*sigmas)/lams
    % b4 = - lams*cos(psi)*(sigmas - sigmat) - (C*r*D*sigmas)    
    A41 = C*diag(D*sigmas);
    A43 = diag(lams.*sin(psi).*(sigmat-sigmas));
    A44 = diag(lams.*cos(psi))+C*diag(r)*D;
    A45 = -diag(lams.*cos(psi));
    A46 = -C*diag(r.*(D*sigmas)./lams);
    b4 = -C*r.*(D*sigmas)+lams.*cos(psi).*(sigmat-sigmas);

    switch params.strainmeasure
    
    case 'pepicelli'

        % A54 = 1
        % A56 = (K*log(lams*lamt))/(lams^2*lamt) - K/(lams^2*lamt) - G/lams^3
        % A57 = G/lamt^3 - K/(lams*lamt^2) + (K*log(lams*lamt))/(lams*lamt^2)
        % b5 = gamma - sigmas - (G*(1/lams^2 - 1/lamt^2))/2 + (K*log(lams*lamt))/(lams*lamt)
    
        % A65 = 1
        % A66 = G/lams^3 - K/(lams^2*lamt) + (K*log(lams*lamt))/(lams^2*lamt)
        % A67 = (K*log(lams*lamt))/(lams*lamt^2) - K/(lams*lamt^2) - G/lamt^3
        % b6 = gamma - sigmat + (G*(1/lams^2 - 1/lamt^2))/2 + (K*log(lams*lamt))/(lams*lamt)

        Asubs = diag((1.-log(lams.*lamt))./(lams.^2));
        Asubt = diag((1.-log(lams.*lamt))./(lamt.^2));
        
        % determine sigma^r
        A55 = -eye(N);
        A56 = Asubs.*diag(Kmod./lamt) - diag(Gmod.*lams.^(-3));
        A57 = Asubt.*diag(Kmod./lams) + diag(Gmod.*lamt.^(-3));
        b5 = -(params.sigma-sigmat+Kmod*log(lams.*lamt)./(lams.*lamt)+...
        0.5*Gmod*(lams.^(-2)-lamt.^(-2)));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = Asubs.*diag(Kmod./lamt) + diag(Gmod*lams.^(-3));
        A67 = Asubt.*diag(Kmod./lams) - diag(Gmod*lamt.^(-3));
        b6 = -(params.sigma-sigmas+Kmod*log(lams.*lamt)./(lams.*lamt)+...
        0.5*Gmod*(lamt.^(-2)-lams.^(-2)));

    end

    % A71 = 1
    % A77 = -rstar
    % b7 = -r + lamt*rstar
    % determine lambda^r
    A71 = eye(N);
    A77 = -diag(itervars.r0);
    b7 = -r+lamt.*itervars.r0;

    % boundary condition dsigmas/ds(0) = 0
    A41(1,:) = ZL;
    A43(1,:) = ZL;
    A44(1,:) = D(1,:);
    A45(1,:) = ZL;
    A46(1,:) = ZL;
    b4(1) = -D(1,:)*sigmas;

    % boundary condition lams(0)=lamt(0)
    % boundary condition lamt(s0) = 1
    A71(1,:) = ZL;
    A77(1,:) = IDL;
    A71(end,:) = ZL;
    A77(end,:) = fliplr(IDL);
    b7(1) =  -lamt(1)+lams(1);
    b7(end) =  1-lamt(end);


    % % combine matrices
    % A = [[A11,   Z, A13,   Z,    Z, A16,   Z,  Z1]; ...
    %      [  Z, A22, A23,   Z,    Z, A26,   Z,  Z1]; ...
    %      [A31, A32, A33, A34,  A35, A36,   Z, A38]; ...
    %      [A41,   Z, A43, A44,  A45, A46,   Z,  Z1]; ...
    %      [  Z,   Z,   Z,   Z,  A55, A56, A57,  Z1]; ...
    %      [  Z,   Z,   Z, A64,    Z, A66, A67,  Z1]; ...
    %      [A71,   Z,   Z,   Z,    Z,   Z, A77,  Z1]; ...
    %      [A81, Z1', A83,  Z1', Z1', A86, Z1',   0]];
    % 
    %b = [b1;b2;b3;b4;b5;b6;b7;b8];

    A = [[A11,   Z, A13,   Z,    Z, A16,   Z,  Z1]; ...
         [  Z, A22, A23,   Z,    Z, A26,   Z,  Z1]; ...
         [A31, A32, A33, A34,  A35, A36,   Z, A38]; ...
         [A41,   Z, A43, A44,  A45, A46,   Z,  Z1]; ...
         [  Z,   Z,   Z,   Z,    I,   Z,   Z,  Z1]; ...
         [  Z,   Z,   Z,   Z,    Z,   I,   Z,  Z1]; ...
         [A71,   Z,   Z,   Z,    Z,   Z, A77,  Z1]; ...
         [A81, Z1', A83,  Z1', Z1', A86, Z1',   0]];
         %[Z1', Z1', Z1', Z1',  Z1', Z1', Z1',   1]];

    b = [b1;b2;b3;b4;0*b5;0*b6;b7;b8];

    %fprintf('    max(b) = %d\n',max(abs([b5;b6])));

end
