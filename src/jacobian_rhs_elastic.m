function [A,b] = jacobian_rhs_simple(params,itervars)
    
    D = params.D;
    w = params.w;
    r = itervars.r;
    z = itervars.z;
    psi = itervars.psi;
    taus = itervars.taus;
    taup = itervars.taup;
    lams = itervars.lams;
    lamp = itervars.lamp;
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
    A11 = C*D;
    A13 = diag(lams.*sin(psi));
    A16 = diag(-cos(psi));
    A18 = D*r;
    b1 = lams.*cos(psi)-C*D*r;

    % determine z from psi (incl lams)
    A22 = C*D;
    A23 = diag(-lams.*cos(psi));
    A26 = diag(-sin(psi));
    A28 = D*z;
    b2 = lams.*sin(psi)-C*D*z;

    % determine psi from laplace law
    A31 = diag(-taup.*sin(psi)./r.^2);
    A32 = diag(lams);
    A33 = C*diag(taus)*D+diag(taup.*cos(psi)./r);
    A34 = C*diag(D*psi);
    A35 = diag(sin(psi)./r);
    A36 = diag(z-p0+taup.*sin(psi)./r);
    A38 = taus.*(D*psi);
    A39 = -lams;
    b3 =  -(C*taus.*(D*psi)+lams.*(z-p0+taup.*sin(psi)./r));

    if params.compresstype == 1 || ii == 1
        % determine pressure - use volume      
        A91 = 2*w.*r'.*sin(psi').*lams';
        A93 = w.*r'.^2.*cos(psi').*lams';
        A98 = -params.volume/pi;
        b9 = -(w*(r.^2.*sin(psi).*lams)-C*params.volume/pi);
    else
        % determine pressure - use area
        A91 = 2*w.*lams';
        A93 = zeros(1,N);
        A98 = -params.area/pi;
        b9 = -(2*w*(r.*lams)-C*params.area/pi);
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

    % determine taus from projection of force balance
    % THIS MUST ME CHECKED: SHOULD EVERY DERIVATIVE D BE MULTIPLIED BY C??!!
    A41 = C*diag(D*taus);
    A43 = diag(lams.*sin(psi).*(taup-taus));
    A44 = diag(lams.*cos(psi))+C*diag(r)*D;
    A46 = diag(cos(psi).*(taup-taus));
    A45 = diag(-lams.*cos(psi));
    b4 = -C*r.*(D*taus)+lams.*cos(psi).*(taup-taus);  % check this eq.

    switch params.strainmeasure
    
    case 'generic'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = diag(Kmod./lams - Gmod.*lams.^(-3));
        A57 = diag(Kmod./lamp + Gmod.*lamp.^(-3));
        b5 = -(params.sigma-taup+Kmod*log(lams.*lamp)+...
        0.5*Gmod*(lams.^(-2)-lamp.^(-2)));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = diag(Kmod./lams + Gmod*lams.^(-3));
        A67 = diag(Kmod./lamp - Gmod*lamp.^(-3));
        b6 = -(params.sigma-taus+Kmod*log(lams.*lamp)+...
        0.5*Gmod*(lamp.^(-2)-lams.^(-2)));
    
    case 'knoche'
    
        % determine sigma^r
        A55 = -diag(lams);
        A56 = diag((Kmod-Gmod)+(params.sigma-taup));
        A57 = (Kmod+Gmod)*eye(N);
        b5 = -((Kmod+Gmod)*(lamp-1)+(Kmod-Gmod)*...
        (lams-1)+lams.*(params.sigma-taup));
        
        % determine lambda^s
        A64 = -diag(lamp);
        A66 = (Kmod+Gmod)*eye(N);
        A67 = diag((Kmod-Gmod)+(params.sigma-taus));
        b6 = -((Kmod+Gmod)*(lams-1)+(Kmod-Gmod)*...
        (lamp-1)+lamp.*(params.sigma-taus));
    
    case 'hookean'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = (Kmod-Gmod)*eye(N);
        A57 = (Kmod+Gmod)*eye(N);
        b5 = -((Kmod+Gmod)*(lamp-1)+(Kmod-Gmod)*...
        (lams-1)+params.sigma-taup);
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = (Kmod+Gmod)*eye(N);
        A67 = (Kmod-Gmod)*eye(N);
        b6 = -((Kmod+Gmod)*(lams-1)+(Kmod-Gmod)*...
        (lamp-1)+params.sigma-taus);
    
    case 'hencky'
    
        % determine sigma^r
        A55 = -eye(N);
        %           A56 = (Kmod-Gmod)*eye(N); % incorrect in code Nagel?
        %           A57 = (Kmod+Gmod)*eye(N);
        A56 = (Kmod-Gmod)*eye(N).*diag(1./lams);
        A57 = (Kmod+Gmod)*eye(N).*diag(1./lamp);          
        b5 = -(Kmod*log(lams.*lamp)+Gmod*log(lamp./lams)+...
        (params.sigma-taup));
        
        % determine lambda^s
        A64 = -eye(N);
        %           A66 = (Kmod+Gmod)*eye(N); % incorrect in code Nagel?
        %           A67 = (Kmod-Gmod)*eye(N);
        A66 = (Kmod+Gmod)*eye(N).*diag(1./lams);
        A67 = (Kmod-Gmod)*eye(N).*diag(1./lamp);          
        b6 = -(Kmod*log(lams.*lamp)+Gmod*log(lams./lamp)+...
        (params.sigma-taus));
    
    case 'pepicelli'
    
        Asubs = diag((1.-log(lams.*lamp))./(lams.^2));
        Asubr = diag((1.-log(lams.*lamp))./(lamp.^2));
        
        % determine sigma^r
        A55 = -eye(N);
        A56 = Asubs.*diag(Kmod./lamp) - diag(Gmod.*lams.^(-3));
        A57 = Asubr.*diag(Kmod./lams) + diag(Gmod.*lamp.^(-3));
        b5 = -(params.sigma-taup+Kmod*log(lams.*lamp)./(lams.*lamp)+...
        0.5*Gmod*(lams.^(-2)-lamp.^(-2)));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = Asubs.*diag(Kmod./lamp) + diag(Gmod*lams.^(-3));
        A67 = Asubr.*diag(Kmod./lams) - diag(Gmod*lamp.^(-3));
        b6 = -(params.sigma-taus+Kmod*log(lams.*lamp)./(lams.*lamp)+...
        0.5*Gmod*(lamp.^(-2)-lams.^(-2)));
        
    case 'balemans'
    
        % determine sigma^r
        A55 = -eye(N);
        A56 = diag(Kmod./lams) - diag(Gmod.*lamp./lams.^(-2));
        A57 = diag(Kmod./lamp) + diag(Gmod./lams);
        b5 = -(params.sigma-taup+Kmod*log(lams.*lamp)+...
        Gmod*(lamp./lams-1));
        
        % determine lambda^s
        A64 = -eye(N);
        A66 = diag(Kmod./lams) + diag(Gmod./lamp);
        A67 = diag(Kmod./lamp) - diag(Gmod*lams./lamp.^(-2));
        b6 = -(params.sigma-taus+Kmod*log(lams.*lamp)+...
        Gmod*(lams./lamp-1));

    end

    % determine lambda^r
    A71 = eye(N);
    A77 = diag(-itervars.r0);
    b7 = -r+lamp.*itervars.r0;
    
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
    
    b4(1) = -D(1,:)*taus;
    b7(1) =  -lamp(1)+lams(1);
    b7(end) =  1-lamp(end);

    % combine matrices
    A = [[A11, Z, A13, Z, Z, A16, Z,Z1]; ...
       [Z, A22, A23, Z, Z, A26, Z, Z1];...
       [A31, A32, A33, A34, A35, A36, Z, A39]; ...
       [A41, Z, A43, A44, A45, A46, Z, Z1]; ...
       [Z, Z, Z, Z, A55, A56, A57, Z1]; ...
       [Z, Z, Z, A64, Z, A66, A67, Z1]; ...
       [A71, Z, Z, Z, Z, Z, A77, Z1]; ...
       [A91, Z1',A93,Z1',Z1',Z1',Z1',0]];
    
    b = [b1;b2;b3;b4;b5;b6;b7;b9];


end