function [ c , f ] = fchebt( y, max_order, shift_order )
% [coefficient, filtered function] = fchebt(function, max_order, shift_order)
% is a Fast Chebychev Transform, based on first kind Chebyshev Function T.
% The function `y` is approximated up to `max_order`. When `max_order` is
% negative the order is optimized based on the performance of each new
% order. `shift_order` will push the optimal order.


% Chebyshev polynomials have singular end points but are better weighted
% (modes not increasing with wavenumber)
% pass 3 times and make a final error check showed better precision

global g_error
% if maximum order is negative the response is computed as f

f = 0;
N = length(y);

% Matlab Chebyshev is too slow
YC = fcbt('cheb',N,abs(max_order)+1);

if size(y,1)==1
    y = y';
end;

% backup function y for final comparison
ys0 = y;

% initialize arguments if not passed
if ~exist('max_order','var')
    max_order = N/2;
end

if ~exist('shift_order','var')
    shift_order = 0;
end

% revarg detemines if the truncation error is evaluated
revarg = (max_order<0);
max_order = abs(max_order);
ds = 2/(N-1);      % element width

% create differentiation, integration and ordinate
[d, ~, w, s] = dif1D('fd',-1+ds,2-2*ds,N-2,5); % numerical integral only on regular strip

% w does overall integration
% create w2 for point wise integration
d(1,:) = [1, zeros(1,N-3)];
w2 = inv(d);
w2 = w2(end,:);
w2(1) = 0;

c = zeros(max_order+1,1);
%tic;
%se = [-1; s; 1];
%we = [w(1:10),w(10:11),w(11:end)];

% 1st pass over 6 modes, c are the modes, from y the fitted part is subtracted
for k = 0:min(6,max_order)
    % 2nd order singular integration analytically
    y0 = YC(1,k+1)*y(1);
    dy = (YC(2,k+1)*y(2)-YC(1,k+1)*y(1))/ds;
    singular_left = (y0+dy)*acos(1-ds)-dy*sqrt(ds*(2-ds));
    y0 = YC(end,k+1)*y(end);
    dy = (YC(end,k+1)*y(end)-YC(end-1,k+1)*y(end-1))/ds;
    singular_right = (y0-dy)*acos(1-ds)+dy*sqrt(ds*(2-ds));
    
    % numerical integral and singular endpoint contribution
    c(k+1) = 2/pi*(w2*(y(2:end-1)./sqrt(1-s.^2).*YC(2:end-1,k+1)) + singular_left+singular_right)/(1+(k==0));
    y = y - c(k+1)*YC(:,k+1);
end

% 2nd pass till 12 modes
for k = 0:min(12,max_order)
    y0 = YC(1,k+1)*y(1);
    dy = (YC(2,k+1)*y(2)-YC(1,k+1)*y(1))/ds;
    singular_left = (y0+dy)*acos(1-ds)-dy*sqrt(ds*(2-ds));
    y0 = YC(end,k+1)*y(end);
    dy = (YC(end,k+1)*y(end)-YC(end-1,k+1)*y(end-1))/ds;
    singular_right = (y0-dy)*acos(1-ds)+dy*sqrt(ds*(2-ds));
    
    tmp = 2/pi*(w*(y(2:end-1)./sqrt(1-s.^2).*YC(2:end-1,k+1)) + singular_left+singular_right)/(1+(k==0));
    c(k+1) = c(k+1)+tmp;
    y = y - tmp*YC(:,k+1);
end

% 3rd pass over all modes
for k = 0:max_order
    y0 = YC(1,k+1)*y(1);
    dy = (YC(2,k+1)*y(2)-YC(1,k+1)*y(1))/ds;
    singular_left = (y0+dy)*acos(1-ds)-dy*sqrt(ds*(2-ds));
    y0 = YC(end,k+1)*y(end);
    dy = (YC(end,k+1)*y(end)-YC(end-1,k+1)*y(end-1))/ds;
    singular_right = (y0-dy)*acos(1-ds)+dy*sqrt(ds*(2-ds));
    
    tmp = 2/pi*(w*(y(2:end-1)./sqrt(1-s.^2).*YC(2:end-1,k+1)) + singular_left+singular_right)/(1+(k==0));
    c(k+1) = c(k+1)+tmp;
    y = y - tmp*YC(:,k+1);
end

c2 = zeros(max_order+1,1);

% 4th pass over all modes (several passes showed lower error)
for k = 0:max_order
    y0 = YC(1,k+1)*y(1);
    dy = (YC(2,k+1)*y(2)-YC(1,k+1)*y(1))/ds;
    singular_left = (y0+dy)*acos(1-ds)-dy*sqrt(ds*(2-ds));
    y0 = YC(end,k+1)*y(end);
    dy = (YC(end,k+1)*y(end)-YC(end-1,k+1)*y(end-1))/ds;
    singular_right = (y0-dy)*acos(1-ds)+dy*sqrt(ds*(2-ds));
    
    c2(k+1) = 2/pi*(w*(y(2:end-1)./sqrt(1-s.^2).*YC(2:end-1,k+1)) + singular_left+singular_right)/(1+(k==0));
end

c = c + c2;

g_error = -1;
errmin = 1e12;
if revarg
    errfind = zeros(max_order+1,1);
    % reconstruct the function as f and determine minimal error
    for k = 0:max_order
        f = f +c(k+1)*YC(:,k+1);
        errfind(k+1) = rms(f-ys0);
        errmin = min(errfind(k+1),errmin);
    end

    % look for the last mode that brought the error within 10% of that
    k = find(errfind<1.5*errmin,1,'first');
    
    % shift order if possible/applicable
    k = min(k+shift_order,length(c)-1);
    
    % erase all higher modes
    c(k+1:end) = 0;
    g_error = errfind;
end

return