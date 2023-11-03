function [d,dd,w,s]=dif1D(type,s0,L,N,pts)
% [d/ds, d/ds^2, integral, s] = dif1D('type',s0,length,N_dof,order), creates
% the 1D differentiation matrices
% These functions have been collected from Jerome Hoepffners teaching
% materials at http://basilisk.fr/sandbox/easystab/README

switch type
 case 'fd' % finite difference 
  scale=L/2;
  [~,d] = fddif(N,1,pts); 
  [s,dd] = fddif(N,2,pts); 
  s=s*scale; 
  d=(d/scale); dd=dd/scale^2; s=s-s(1)+s0;
    d = full(d);
  w=([diff(s'),0]+[0,diff(s')])/2;  
 case 'cheb' % chebychev 
  scale=-L/2;
  [s,DM] = chebdif(N,2);
  d=DM(:,:,1);  
  dd=DM(:,:,2);
  s=s*scale; 
  d=d/scale; dd=dd/scale^2; s=s-s(1)+s0;
  w=L*clencurt(N)/2; 
end

end

function [x,D]=fddif(N,order,pts)
% build equispaced grid on [-1,1], and 
% five points finite diference matrix for N mesh points 

%pts=5;
x=linspace(-1,1,N)';
h=x(2)-x(1);

% subroutine for finite difference weights
W=ufdwt(h,pts,order);
t=(pts+1)/2;

%central difference in the middle
D=spdiags(ones(N,1)*W(t,:),-t+1:t-1,N,N);

for indd=1:t-1
  D(indd,1:pts)=W(indd,:);   % at the left boundary
  D(N-indd+1,end-pts+1:end)=W(end-indd+1,:); % at the right boundary
end

end

% Diffm compute D = differentiation matrix, x = Chebyshev grid
% Order is N, in the interval a to b
%  function [D,x] = dmat(N,a,b)
%  if N==0, D=0; x=1; return, end
%  if a>=b, D=0; x=0; disp('point a smaller than b'); return, end 
%  x = (a-b)/2*cos(pi*(0:N)/N)' + (a+b)/2;
%  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
%  X = repmat(x,1,N+1);
%  dX = X-X';                  
%  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
%  D  = D - diag(sum(D,2));                 % diagonal entries
%  end

function W=ufdwt(h,pts,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ufdwt.m
%
% Compute Finite Difference Weights for a Uniform Grid
%
% Input Parameters:
%
% h     - spacing between FD nodes
% pts   - number of FD points in scheme (3-pt, 5-pt, etc)
% order - order of the derivative operator to compute weights for
%         (note: order<pts-1!)
%         1 computes first derivative differences       
%         2 computes second derivative differences, etc
%  
% Output Parameter:
%
% W is the weight matrix. Each row contains a different set of weights
% (centered or off). If, for example, the number of finite difference points
% is odd, the centered difference weights will appear in the middle row.
%
% Written by: Greg von Winckel - 06/16/04
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=2*pts-1;  p1=pts-1;

A=repmat((0:p1)',1,N);      
B=repmat((-p1:p1)*h,pts,1);

M=(B.^A)./gamma(A+1); 

rhs=zeros(pts,1);   rhs(order+1)=1;

W=zeros(pts,pts);

for k=1:pts
    W(:,k)=M(:,(0:p1)+k)\rhs;
end

W=W';   W(1:pts,:)=W(pts:-1:1,:);

end

function DM = poldif(x, malpha, B)

%  The function DM =  poldif(x, maplha, B) computes the
%  differentiation matrices D1, D2, ..., DM on arbitrary nodes.
%
%  The function is called with either two or three input arguments.
%  If two input arguments are supplied, the weight function is assumed 
%  to be constant.   If three arguments are supplied, the weights should 
%  be defined as the second and third arguments.
%
%  Input (constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   M, the number of derivatives required (integer).
%  B:        Omitted.
%
%  Note:     0 < M < N-1.
%
%  Input (non-constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
%  B:        Matrix of size M x N,  where M is the highest 
%            derivative required.  It should contain the quantities 
%            B(ell,j) = beta(ell,j) = (ell-th derivative
%            of alpha(x))/alpha(x),   evaluated at x = x(j).
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.

%  J.A.C. Weideman, S.C. Reddy 1998

       N = length(x);                      
       x = x(:);                     % Make sure x is a column vector 

if nargin == 2                       % Check if constant weight function
       M = malpha;                   % is to be assumed.
   alpha = ones(N,1);              
       B = zeros(M,N);
elseif nargin == 3
   alpha = malpha(:);                % Make sure alpha is a column vector
       M = length(B(:,1));           % First dimension of B is the number 
end                                  % of derivative matrices to be computed
 
        I = eye(N);                  % Identity matrix.
        L = logical(I);              % Logical identity matrix.

       XX = x(:,ones(1,N));
       DX = XX-XX';                  % DX contains entries x(k)-x(j).

    DX(L) = ones(N,1);               % Put 1's one the main diagonal.

        c = alpha.*prod(DX,2);       % Quantities c(j).

        C = c(:,ones(1,N)); 
        C = C./C';                   % Matrix with entries c(k)/c(j).
   
        Z = 1./DX;                   % Z contains entries 1/(x(k)-x(j))
     Z(L) = zeros(N,1);              % with zeros on the diagonal.

        X = Z';                      % X is same as Z', but with 
     X(L) = [];                      % diagonal entries removed.
        X = reshape(X,N-1,N);

        Y = ones(N-1,N);             % Initialize Y and D matrices.
        D = eye(N);                  % Y is matrix of cumulative sums,
        DM = zeros(N,N,M);                             % D differentiation matrices.
for ell = 1:M
        Y   = cumsum([B(ell,:); ell*Y(1:N-1,:).*X]); % Diagonals
        D   = ell*Z.*(C.*repmat(diag(D),1,N) - D);   % Off-diagonals
        D(L)   = Y(N,:);                                % Correct the diagonal
        DM(:,:,ell) = D;                                     % Store the current D
end

end

function IW=clencurt(N)
%
% Computes the integration weigths for pseudo-chebychev on domain [-1 1]
%
% INPUTS:
% N  : the number of points 
%
% OUTPUT:
% IW : vector of the integration weigths.

nW=0:1:N-1;
jW=0:1:N-1;

bW=ones(1,N); bW(1)=0.5; bW(N)=0.5;
cW=2*bW;
bW=bW/(N-1);

S=cos(nW(3:N)'*jW*(pi/(N-1)));
IW=bW.*( (2+(cW(3:N).*((1+(-1).^nW(3:N))./(1-nW(3:N).^2)))*S) );

end

function [x, DM] = chebdif(N, M)

%  The function DM =  chebdif(N,M) computes the differentiation 
%  matrices D1, D2, ..., DM on Chebyshev nodes. 
% 
%  Input:
%  N:        Size of differentiation matrix.        
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
    
%  J.A.C. Weideman, S.C. Reddy 1998.

     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

    n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

     k = (0:N-1)';                        % Compute theta vector.
    th = k*pi/(N-1);

     x = sin(pi*(N-1:-2:1-N)'/(2*(N-1))); % Compute Chebyshev points.

     T = repmat(th/2,1,N);                
    DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity. 
    DX = [DX(1:n1,:); -rot90(DX(1:n2,:),2)];   % Flipping trick. 
 DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.

     C = toeplitz((-1).^k);               % C is the matrix with 
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
  Z(L) = zeros(N,1);                      % with zeros on the diagonal.

     D = eye(N);                          % D contains diff. matrices.
     DM = zeros(N,N,M);                                     
for ell = 1:M
          D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
       D(L) = -sum(D,2)';                            % Correct main diagonal of D
DM(:,:,ell) = D;                                   % Store current D in DM
end
end
