% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area
% NOTE: the input is dimensionfull, but the computation is dimensionless
% after the computation, the variables need to be properly scaled back!

close all; clear

addpath('subs/')

% physical parameters
sigma = 72;        % surface tension [mN/m]
grav = 9.807e3;    % gravitational acceleration [mm/s^2]
rneedle = 1;       % radius of the needle [mm]
volume0 = 32;      % prescribed volume in mm^3
deltarho = 1e-3;   % density difference [10^6 kg/m^3]

% numerical parameters
N = 40;            % resolution of the discretization for calculation
Nplot = 80;        % resolution of the discretization for plotting
Ncheb = 10;        % number of Chebyshev to describe the shape
alpha = 1;       % relaxation parameter in the Newton-Raphson scheme
Kmod = 20;         % elastic dilational modulus [mN/m]
Gmod = 20;         % elastic shear modulus [mN/m]
compresstype = 1;  % 1: compress the volume    other: compress the area
frac = [0.6];      % compute elastic stresses for these compressions
strainmeasure = 'pepicelli'; % which elastic constitutive model

% calculate the dimensionless quantities
sigmaprime = sigma/(deltarho*grav*rneedle^2);
volume0prime = volume0/rneedle^3;
Kmodprime = Kmod/(deltarho*grav*rneedle^2);
Gmodprime = Gmod/(deltarho*grav*rneedle^2);

% open figure for plotting the shape
figure; hold on;
rmax = 0; zmin = 1e12; % variables used for plotting

for ii = 1:length(frac)+1
    
  % first iteration is done isotropic
  elastic = (ii ~= 1);  % elastic = 0: simple Young-Laplace interface       
                        % elastic = 1: elastic interface

  % determine the current target volume
  if ~elastic
    % predict the maximum length of the interface (empirical Nagel)
    smax = sqrt(sigmaprime)*2.0/0.8701;
    volumeprime = volume0prime;
  else
    % volume/area given by compression ratios
    if compresstype == 1
      volume = volume0*frac(ii-1);
      volumeprime = volume/rneedle^3;
    else
      areaprime = areaprevious*frac(ii-1);
    end
  end

  % get the differentation/integration matrices and the grid
  [D,~,w,s] = dif1D('cheb',0,smax,N,5);

  % predict the shape of the interface (empirical Nagel)
  z = -4/3*smax/pi*(cos(pi*3/4*s/smax));
  z = z - max(z);
  r = 4/3*smax/pi*(sin(pi*3/4*s/smax));
  psi = pi*3/4*s/smax;

  % initialize the surface strains ans tresses
  lamp = ones(N,1); lams = lamp;
  taus = sigmaprime*ones(N,1); taup = taus;

  if ~elastic
    C = 1; % initial stretch parameter
    r0 = r;
  end

  p0 = sqrt(sigmaprime)*1.5; % predict the pressure (empirical Nagel)

  % initialize some variables 
  Z = zeros(N);            % matrix filled with zeros
  IDL = [1, zeros(1,N-1)]; % line with single one and rest zeros
  ZL = zeros(1,N);         % line completely filled with zeros
  u = ones(3*N+2,1); b = ones(3*N+2,1); % solution vector and RHS
  iter = 0; crash = 0; 

  while rms(u) > 1e-10

    iter = iter + 1;
    
    if iter > 1200 
      error('iter > 1200!');
    end

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
    
    if compresstype == 1 || ii == 1
      % determine pressure - use volume      
      A91 = 2*w.*r'.*sin(psi').*lams';
      A93 = w.*r'.^2.*cos(psi').*lams';
      A98 = -volumeprime/pi;
      b9 = -(w*(r.^2.*sin(psi).*lams)-C*volumeprime/pi);
    else
      % determine pressure - use area
      A91 = 2*w.*lams';
      A93 = zeros(1,N);
      A98 = -areaprime/pi;
      b9 = -(2*w*(r.*lams)-C*areaprime/pi);
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

    if ~elastic

      % determine arclength scale
      A81 = fliplr(IDL);
      b8 = (1-r(end));

      A = [[A11, Z, A13, A18, Z1];[Z, A22, A23, A28, Z1];
          [A31, A32, A33, A38, A39];[A81, zeros(1,2*N), 0,0];
          [A91, Z1',A93,A98,0]];
      b = [b1;b2;b3;b8;b9];

    else

      % determine taus from projection of force balance
      % THIS MUST ME CHECKED: SHOULD EVERY DERIVATIVE D BE MULTIPLIED BY C??!!
      A41 = C*diag(D*taus);
      A43 = diag(lams.*sin(psi).*(taup-taus));
      A44 = diag(lams.*cos(psi))+C*diag(r)*D;
      A46 = diag(cos(psi).*(taup-taus));
      A45 = diag(-lams.*cos(psi));
      b4 = -C*r.*(D*taus)+lams.*cos(psi).*(taup-taus);  % check this eq.

      switch strainmeasure

        case 'generic'

          % determine sigma^r
          A55 = -eye(N);
          A56 = diag(Kmodprime./lams - Gmodprime.*lams.^(-3));
          A57 = diag(Kmodprime./lamp + Gmodprime.*lamp.^(-3));
          b5 = -(sigmaprime-taup+Kmodprime*log(lams.*lamp)+...
            0.5*Gmodprime*(lams.^(-2)-lamp.^(-2)));

          % determine lambda^s
          A64 = -eye(N);
          A66 = diag(Kmodprime./lams + Gmodprime*lams.^(-3));
          A67 = diag(Kmodprime./lamp - Gmodprime*lamp.^(-3));
          b6 = -(sigmaprime-taus+Kmodprime*log(lams.*lamp)+...
            0.5*Gmodprime*(lamp.^(-2)-lams.^(-2)));

        case 'knoche'

          % determine sigma^r
          A55 = -diag(lams);
          A56 = diag((Kmodprime-Gmodprime)+(sigmaprime-taup));
          A57 = (Kmodprime+Gmodprime)*eye(N);
          b5 = -((Kmodprime+Gmodprime)*(lamp-1)+(Kmodprime-Gmodprime)*...
            (lams-1)+lams.*(sigmaprime-taup));

          % determine lambda^s
          A64 = -diag(lamp);
          A66 = (Kmodprime+Gmodprime)*eye(N);
          A67 = diag((Kmodprime-Gmodprime)+(sigmaprime-taus));
          b6 = -((Kmodprime+Gmodprime)*(lams-1)+(Kmodprime-Gmodprime)*...
            (lamp-1)+lamp.*(sigmaprime-taus));

        case 'hookean'

          % determine sigma^r
          A55 = -eye(N);
          A56 = (Kmodprime-Gmodprime)*eye(N);
          A57 = (Kmodprime+Gmodprime)*eye(N);
          b5 = -((Kmodprime+Gmodprime)*(lamp-1)+(Kmodprime-Gmodprime)*...
            (lams-1)+sigmaprime-taup);

          % determine lambda^s
          A64 = -eye(N);
          A66 = (Kmodprime+Gmodprime)*eye(N);
          A67 = (Kmodprime-Gmodprime)*eye(N);
          b6 = -((Kmodprime+Gmodprime)*(lams-1)+(Kmodprime-Gmodprime)*...
            (lamp-1)+sigmaprime-taus);

        case 'hencky'

          % determine sigma^r
          A55 = -eye(N);
%           A56 = (Kmodprime-Gmodprime)*eye(N); % incorrect in code Nagel?
%           A57 = (Kmodprime+Gmodprime)*eye(N);
          A56 = (Kmodprime-Gmodprime)*eye(N).*diag(1./lams);
          A57 = (Kmodprime+Gmodprime)*eye(N).*diag(1./lamp);          
          b5 = -(Kmodprime*log(lams.*lamp)+Gmodprime*log(lamp./lams)+...
            (sigmaprime-taup));

          % determine lambda^s
          A64 = -eye(N);
%           A66 = (Kmodprime+Gmodprime)*eye(N); % incorrect in code Nagel?
%           A67 = (Kmodprime-Gmodprime)*eye(N);
          A66 = (Kmodprime+Gmodprime)*eye(N).*diag(1./lams);
          A67 = (Kmodprime-Gmodprime)*eye(N).*diag(1./lamp);          
          b6 = -(Kmodprime*log(lams.*lamp)+Gmodprime*log(lams./lamp)+...
            (sigmaprime-taus));

        case 'pepicelli'
            
          Asubs = diag((1.-log(lams.*lamp))./(lams.^2));
          Asubr = diag((1.-log(lams.*lamp))./(lamp.^2));

          % determine sigma^r
          A55 = -eye(N);
          A56 = Asubs.*diag(Kmodprime./lamp) - diag(Gmodprime.*lams.^(-3));
          A57 = Asubr.*diag(Kmodprime./lams) + diag(Gmodprime.*lamp.^(-3));
          b5 = -(sigmaprime-taup+Kmodprime*log(lams.*lamp)./(lams.*lamp)+...
            0.5*Gmodprime*(lams.^(-2)-lamp.^(-2)));

          % determine lambda^s
          A64 = -eye(N);
          A66 = Asubs.*diag(Kmodprime./lamp) + diag(Gmodprime*lams.^(-3));
          A67 = Asubr.*diag(Kmodprime./lams) - diag(Gmodprime*lamp.^(-3));
          b6 = -(sigmaprime-taus+Kmodprime*log(lams.*lamp)./(lams.*lamp)+...
            0.5*Gmodprime*(lamp.^(-2)-lams.^(-2)));
        
        case 'balemans'

          % determine sigma^r
          A55 = -eye(N);
          A56 = diag(Kmodprime./lams) - diag(Gmodprime.*lamp./lams.^(-2));
          A57 = diag(Kmodprime./lamp) + diag(Gmodprime./lams);
          b5 = -(sigmaprime-taup+Kmodprime*log(lams.*lamp)+...
            Gmodprime*(lamp./lams-1));

          % determine lambda^s
          A64 = -eye(N);
          A66 = diag(Kmodprime./lams) + diag(Gmodprime./lamp);
          A67 = diag(Kmodprime./lamp) - diag(Gmodprime*lams./lamp.^(-2));
          b6 = -(sigmaprime-taus+Kmodprime*log(lams.*lamp)+...
            Gmodprime*(lams./lamp-1));
        
      end

      % determine lambda^r
      A71 = eye(N);
      A77 = diag(-r0);
      b7 = -r+lamp.*r0;

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

    % if G is 0 the equations are singular, better use a very small G
    if (Gmodprime~=0)
      u = A\b;
    else
      warning('G is zero (singular equations). Consider increasing value');
      u = A\b;
    end
    
    clear A*

    % update variables
    r   = r + alpha*u(1:N);
    z   = z + alpha*u(N+1:2*N);
    psi = psi + alpha*u(2*N+1:3*N);
    p0  = p0 + alpha*u(end);

    if elastic
      taus = taus + alpha*u(3*N+1:4*N);
      taup = taup + alpha*u(4*N+1:5*N);
      lams = lams + alpha*u(5*N+1:6*N);
      lamp = lamp + alpha*u(6*N+1:7*N);
    else
      C = C+alpha*u(3*N+1);
    end

    if rms(b)>1e3
      crash = 1; break;
    end

    fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

  end

  % compute volume and area
  Vol = pi*w*(r.^2.*sin(psi).*lams)/C;
  Ar = pi*2*w*(r.*lams)/C;
  
  % if isotropic save radius^init, and a^init
  if ~elastic
    r0 = r; z0 = z; areaprevious = Ar;
  end
    
  disp(['volume = ', num2str(rneedle^3*Vol),' mm^3']);
  disp(['area = ', num2str(rneedle^2*Ar),' mm^2']);
  disp(['pressure = ', num2str(deltarho*grav*rneedle*p0),' Pa']);

  % interpolate the numerical solutions on a finer grid. 
  % NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
  % though all the points (see book of Trefethen on Spectral Methods in 
  % Matlab, page  63). For plotting purposes we use a simpler interpolation 
  ss = linspace(s(1),s(end),Nplot)';
  rr = interp1(s,r,ss,'pchip');
  zz = interp1(s,z,ss,'pchip');

  % plot the droplet shape
  plot(rneedle*rr,rneedle*zz); 
  rmax = max([rmax,rr']);
  zmin = min([zmin,zz']);
  
  % rescale the plot
  xlim([0 1.2*rneedle*rmax]);
  ylim([1.2*rneedle*zmin 0]);
  set(gca,'DataAspectRatio',[1 1 1])

  % compute the curvatures (NOTE: d/ds operator is given by C*D, see Nagel)
  kappas = (C*D*psi)./lams;
  kappap = sin(psi)./r;
  kappap(1) = kappas(1);

  % NOTE: there are three coordinates involved: s0 (guessed domain length), 
  % s* (domain for isotropic solution, which is also the reference domain 
  % for the elastic problem) and s (domain in deformed state). The
  % grid (and thus the differentation/integration operators) is defined for 
  % the guessed domain. To obtain the other domains, we
  % use: s* = s0 / C  and  s = \int lambdas ds* = \int lambdas ds / C
  
  % construct the integration matrix from the integration vector
  wmat = repmat(w,N,1);
  wmat = tril(wmat);
  
  % compute the value of s in the deformed state
  sdef = wmat*lams/C;
  
  if crash
    disp(['Calculation crashed, check that initial volume or',...
    'area fraction is okay, then check initial arc-length smax']);
  end
  
end

% calculate the Chebyshev coefficients
coefr = fchebt(r,Ncheb,0);
coefz = fchebt(z,Ncheb,0);
