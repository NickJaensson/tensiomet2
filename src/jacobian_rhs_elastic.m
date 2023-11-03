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

    % Eq. 1-4 Knoche, p85, eq.5.7, Eq. 5-7, eq.5.8
    % determine r from psi (incl lams)
    % A11 = C/lams*D
    % A13 = sin(psi)
    % A16 = -(C*D*r)/lams^2
    % b1 = cos(psi) - (C*D*r)/lams
    A11 = C*D;
    A13 = diag(lams.*sin(psi));
    A16 = diag(-cos(psi));
    A18 = D*r;
    b1 = lams.*cos(psi)-C*D*r;

    % determine z from psi (incl lams)
    % A11 = C/lams*D
    % A13 = sin(psi)
    % A16 = -(C*D*r)/lams^2
    % b1 = cos(psi) - (C*D*r)/lams
    A22 = C*D;
    A23 = diag(-lams.*cos(psi));
    A26 = diag(-sin(psi));
    A28 = D*z;
    b2 = lams.*sin(psi)-C*D*z;

    % determine psi from laplace law
    % A31 = -(sigmat*sin(psi))/r^2
    % A32 = g*rho
    % A33 = (sigmat*cos(psi))/r + (C*sigmas)/lams*D
    % A34 = (C*D*psi)/lams
    % A35 = sin(psi)/r
    % A36 = -(C*D*psi*sigmas)/lams^2
    % A38 = -1
    % b3 = P - (sigmat*sin(psi))/r - g*rho*z - (C*D*psi*sigmas)/lams
    A31 = diag(-sigmat.*sin(psi)./r.^2);
    A32 = diag(lams);
    A33 = C*diag(sigmas)*D+diag(sigmat.*cos(psi)./r);
    A34 = C*diag(D*psi);
    A35 = diag(sin(psi)./r);
    A36 = diag(z-p0+sigmat.*sin(psi)./r);
    A38 = sigmas.*(D*psi);
    A39 = -lams;
    b3 =  -(C*sigmas.*(D*psi)+lams.*(z-p0+sigmat.*sin(psi)./r));

    % A81 = 2*int*r*pi*sin(psi)
    % A83 = int*r^2*pi*cos(psi)
    % b8 = C*V - int*r^2*pi*sin(psi)
    if params.compresstype == 1 || ii == 1
        % determine pressure - use volume      
        A81 = 2*w.*r'.*sin(psi').*lams';
        A83 = w.*r'.^2.*cos(psi').*lams';
        A88 = -params.volume/pi;
        b8 = -(w*(r.^2.*sin(psi).*lams)-C*params.volume/pi);
    else
        % determine pressure - use area
        A81 = 2*w.*lams';
        A83 = zeros(1,N);
        A88 = -params.area/pi;
        b8 = -(2*w*(r.*lams)-C*params.area/pi);
    end

    % Boundary conditions
    A11(1,:) = IDL;
    A13(1,:) = ZL;
    A16(1,:) = ZL;
    A18(1) = 0;
    A22(1,:) = fliplr(IDL);
    A23(1,:) = ZL;
    A26(1,:) = ZL;
    A28(1) = 0;
    A31(1,:) = ZL;
    A32(1,:) = ZL;
    A33(1,:) = IDL;
    A34(1,:) = ZL;
    A35(1,:) = ZL;
    A36(1,:) = ZL;
    A38(1,:) = 0;
    A39(1,:) = 0;

    b1(1) = -r(1);
    b2(1) = -z(end);
    b3(1) = -psi(1);

    Z1 = zeros(N,1);

    % determine sigmas from projection of force balance
    % A41 = (C*D*sigmas)/lams
    % A43 = -sin(psi)*(sigmas - sigmat)
    % A44 = cos(psi) + (C*r)/lams*D
    % A45 = -cos(psi)
    % A46 = -(C*r*D*sigmas)/lams^2
    % b4 = - cos(psi)*(sigmas - sigmat) - (C*r*D*sigmas)/lams    
    A41 = C*diag(D*sigmas);
    A43 = diag(lams.*sin(psi).*(sigmat-sigmas));
    A44 = diag(lams.*cos(psi))+C*diag(r)*D;
    A46 = diag(cos(psi).*(sigmat-sigmas));
    A45 = diag(-lams.*cos(psi));
    b4 = -C*r.*(D*sigmas)+lams.*cos(psi).*(sigmat-sigmas);  % check this eq.

    switch params.strainmeasure
    
    case 'generic'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = diag(Kmod./lams - Gmod.*lams.^(-3));
        A57 = diag(Kmod./lamt + Gmod.*lamt.^(-3));
        b5 = -(params.sigma-sigmat+Kmod*log(lams.*lamt)+...
        0.5*Gmod*(lams.^(-2)-lamt.^(-2)));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = diag(Kmod./lams + Gmod*lams.^(-3));
        A67 = diag(Kmod./lamt - Gmod*lamt.^(-3));
        b6 = -(params.sigma-sigmas+Kmod*log(lams.*lamt)+...
        0.5*Gmod*(lamt.^(-2)-lams.^(-2)));
    
    case 'knoche'
    
        % determine sigma^r
        A55 = -diag(lams);
        A56 = diag((Kmod-Gmod)+(params.sigma-sigmat));
        A57 = (Kmod+Gmod)*eye(N);
        b5 = -((Kmod+Gmod)*(lamt-1)+(Kmod-Gmod)*...
        (lams-1)+lams.*(params.sigma-sigmat));
        
        % determine lambda^s
        A64 = -diag(lamt);
        A66 = (Kmod+Gmod)*eye(N);
        A67 = diag((Kmod-Gmod)+(params.sigma-sigmas));
        b6 = -((Kmod+Gmod)*(lams-1)+(Kmod-Gmod)*...
        (lamt-1)+lamt.*(params.sigma-sigmas));
    
    case 'hookean'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = (Kmod-Gmod)*eye(N);
        A57 = (Kmod+Gmod)*eye(N);
        b5 = -((Kmod+Gmod)*(lamt-1)+(Kmod-Gmod)*...
        (lams-1)+params.sigma-sigmat);
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = (Kmod+Gmod)*eye(N);
        A67 = (Kmod-Gmod)*eye(N);
        b6 = -((Kmod+Gmod)*(lams-1)+(Kmod-Gmod)*...
        (lamt-1)+params.sigma-sigmas);
    
    case 'hencky'
    
        % determine sigma^r
        A55 = -eye(N);
        %           A56 = (Kmod-Gmod)*eye(N); % incorrect in code Nagel?
        %           A57 = (Kmod+Gmod)*eye(N);
        A56 = (Kmod-Gmod)*eye(N).*diag(1./lams);
        A57 = (Kmod+Gmod)*eye(N).*diag(1./lamt);          
        b5 = -(Kmod*log(lams.*lamt)+Gmod*log(lamt./lams)+...
        (params.sigma-sigmat));
        
        % determine lambda^s
        A64 = -eye(N);
        %           A66 = (Kmod+Gmod)*eye(N); % incorrect in code Nagel?
        %           A67 = (Kmod-Gmod)*eye(N);
        A66 = (Kmod+Gmod)*eye(N).*diag(1./lams);
        A67 = (Kmod-Gmod)*eye(N).*diag(1./lamt);          
        b6 = -(Kmod*log(lams.*lamt)+Gmod*log(lams./lamt)+...
        (params.sigma-sigmas));
    
    case 'pepicelli'


        % A54 = 1
        % A56 = (K*log(lams*lamt))/(lams^2*lamt) - K/(lams^2*lamt) - G/lams^3 + 0*D
        % A57 = G/lamt^3 - K/(lams*lamt^2) + (K*log(lams*lamt))/(lams*lamt^2) + 0*D
        % b5 = gamma - sigmas - (G*(1/lams^2 - 1/lamt^2))/2 + (K*log(lams*lamt))/(lams*lamt)
    
        % A65 = 1
        % A66 = G/lams^3 - K/(lams^2*lamt) + (K*log(lams*lamt))/(lams^2*lamt) + 0*D
        % A67 = (K*log(lams*lamt))/(lams*lamt^2) - K/(lams*lamt^2) - G/lamt^3 + 0*D
        % b6 = gamma - sigmat + (G*(1/lams^2 - 1/lamt^2))/2 + (K*log(lams*lamt))/(lams*lamt)

        Asubs = diag((1.-log(lams.*lamt))./(lams.^2));
        Asubr = diag((1.-log(lams.*lamt))./(lamt.^2));
        
        % determine sigma^r
        A55 = -eye(N);
        A56 = Asubs.*diag(Kmod./lamt) - diag(Gmod.*lams.^(-3));
        A57 = Asubr.*diag(Kmod./lams) + diag(Gmod.*lamt.^(-3));
        b5 = -(params.sigma-sigmat+Kmod*log(lams.*lamt)./(lams.*lamt)+...
        0.5*Gmod*(lams.^(-2)-lamt.^(-2)));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = Asubs.*diag(Kmod./lamt) + diag(Gmod*lams.^(-3));
        A67 = Asubr.*diag(Kmod./lams) - diag(Gmod*lamt.^(-3));
        b6 = -(params.sigma-sigmas+Kmod*log(lams.*lamt)./(lams.*lamt)+...
        0.5*Gmod*(lamt.^(-2)-lams.^(-2)));
        
    case 'balemans'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = diag(Kmod./lams) - diag(Gmod.*lamt./lams.^(-2));
        A57 = diag(Kmod./lamt) + diag(Gmod./lams);
        b5 = -(params.sigma-sigmat+Kmod*log(lams.*lamt)+...
        Gmod*(lamt./lams-1));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = diag(Kmod./lams) + diag(Gmod./lamt);
        A67 = diag(Kmod./lamt) - diag(Gmod*lams./lamt.^(-2));
        b6 = -(params.sigma-sigmas+Kmod*log(lams.*lamt)+...
        Gmod*(lams./lamt-1));

    end

    % A65 = 1
    % A66 = G/lams^3 - K/(lams^2*lamt) + (K*log(lams*lamt))/(lams^2*lamt) + 0*D
    % A67 = (K*log(lams*lamt))/(lams*lamt^2) - K/(lams*lamt^2) - G/lamt^3 + 0*D
    % b6 = gamma - sigmat + (G*(1/lams^2 - 1/lamt^2))/2 + (K*log(lams*lamt))/(lams*lamt)
    % determine lambda^r
    A71 = eye(N);
    A77 = diag(-itervars.r0);
    b7 = -r+lamt.*itervars.r0;
    
    % Boundary conditions
    A41(1,:) = ZL;
    A43(1,:) = ZL;
    A44(1,:) = D(1,:);
    A45(1,:) = ZL;
    A46(1,:) = ZL;
    
    A71(1,:) = ZL;
    A77(1,:) = IDL;
    A71(end,:) = ZL;
    A77(end,:) = fliplr(IDL);
    
    b4(1) = -D(1,:)*sigmas;
    b7(1) =  -lamt(1)+lams(1);
    b7(end) =  1-lamt(end);

    % combine matrices
    A = [[A11, Z, A13, Z, Z, A16, Z,Z1]; ...
       [Z, A22, A23, Z, Z, A26, Z, Z1];...
       [A31, A32, A33, A34, A35, A36, Z, A39]; ...
       [A41, Z, A43, A44, A45, A46, Z, Z1]; ...
       [Z, Z, Z, Z, A55, A56, A57, Z1]; ...
       [Z, Z, Z, A64, Z, A66, A67, Z1]; ...
       [A71, Z, Z, Z, Z, Z, A77, Z1]; ...
       [A81, Z1',A83,Z1',Z1',Z1',Z1',0]];
    
    b = [b1;b2;b3;b4;b5;b6;b7;b8];


end