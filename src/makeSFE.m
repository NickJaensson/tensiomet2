function [GK, lams, lamr] = makeSFE(strainmeasure,vars_sol_ref,vars_num_ref,vars_sol,vars_num,params_num,verbose)
% makeSFE(guessed G and K) = [[GK], lams, lamr]
% Computes the deformations with a material model and returns the material
% parameters. Does not need an initial state or a surfac tension

% g_strainmeasure: strain measure, e.g., hencky
% g_memptr: g_memptr(state): pointer to state=state1 and state=state2
% g_error: error of the fitted moduli
% g_echo: print some info to screen or not? 
% glob_w(state,:): w in state=1 and state=2 (determined in getShape)
% glob_d(:,:,state): d in state=1 and state=2 (determined in getShape)
% glob_s(:,state): s in state=1 and state=2 (determined in getShape)
% glob_r(:,state): r in state=1 and state=2 (determined in getShape)
% glob_z(:,state): z in state=1 and state=2 (determined in getShape)
% glob_ts(:,state): taus in state=1 and state=2 (determined in makeCMD)
% glob_tr(:,state): taup in state=1 and state=2 (determined in makeCMD)

g_strainmeasure = strainmeasure;
g_memptr = [1,2];

glob_w = zeros(2,params_num.N);
glob_d = zeros(params_num.N,params_num.N,2);
glob_s = zeros(params_num.N,2);
glob_r = zeros(params_num.N,2);
glob_z = zeros(params_num.N,2);
glob_ts = zeros(params_num.N,2);
glob_tr = zeros(params_num.N,2);

% BELOW SHOULD ALL BE CHECKED TO MAKE SURE WE USE THE CORRECT NUMERICAL
% VARIABLES !

glob_w(1,:) = vars_num_ref.ws;
glob_d(:,:,1) = vars_num_ref.Ds;
glob_s(:,1) = vars_num_ref.s;
glob_r(:,1) = vars_sol_ref.r;
glob_z(:,1) = vars_sol_ref.z;
glob_ts(:,1) = vars_sol_ref.sigmas;
glob_tr(:,1) = vars_sol_ref.sigmap;

glob_w(2,:) = vars_num.ws;
glob_d(:,:,2) = vars_num.Ds;
glob_s(:,2) = vars_num.s;
glob_r(:,2) = vars_sol.r;
glob_z(:,2) = vars_sol.z;
glob_ts(:,2) = vars_sol.sigmas;
glob_tr(:,2) = vars_sol.sigmap;


if length(g_memptr)<2
    disp('Requires g_memptr to indicate 2 states');
    GK= [-1,-1];
    lams = -ones(length(glob_r));
    lamr = -ones(length(glob_r)); 
    return
end

ptr1 = g_memptr(1);
ptr2 = g_memptr(2);

% initialize variables
wold = glob_w(ptr2,:);
dold = glob_d;
sold = glob_s;
rold = glob_r;
zold = glob_z;
ts = glob_ts;
tr = glob_tr;

N = length(tr(:,1));

K = params_num.K_guess;
G = params_num.G_guess;

% reference state variables and derivatives
dsds0 = sold(end,ptr2)/sold(end,ptr1);
lams = ones(N,1)*dsds0;
lamr = rold(:,ptr2)./rold(:,ptr1);
lamr(1) = lams(1);
u = 1e12*[1;1];

iter = 1;
oldb = 1e12;
d = dold(:,:,ptr2);
di = d;
di(1,:) = [1,zeros(1,N-1)];
Itg = inv(di);
Itg(:,1) = zeros(N,1);

% reference state stresses for interpolation
s1 = ts(:,ptr1);
s2 = tr(:,ptr1);

u=ones(N+2,1); % initialize to get while-loop going

while rms(u(N+1:N+2))>params_num.eps_inv_sfe


    if iter > params_num.maxiter_inv_sfe
        error('Iteration did not converge!')
    end  

    % material equations
    if strcmp(g_strainmeasure,'linear_hookean')

        A  = [[            eye(N)*2*K, 2*(lams+lamr-2),  zeros(N,1)];...
              [            eye(N)*2*G,      zeros(N,1), 2*(lams-lamr)];...
              [-wold*diag(1./lams.^2),               0,          0]];

        b = [ts(:,ptr2)+tr(:,ptr2)-ts(:,ptr1)-tr(:,ptr1)-2*K*(lams+lamr-2);...
                               ts(:,ptr2)-tr(:,ptr2)-s1+s2-2*G*(lams-lamr);...
                                             sold(end,ptr1)-wold*(1./lams)];
    elseif strcmp(g_strainmeasure,'hencky')

        A  = [[       diag(2*K./lams), 2*log(lams.*lamr),        zeros(N,1)];...
              [       diag(2*G./lams),        zeros(N,1), 2*log(lams./lamr)];...
              [-wold*diag(1./lams.^2),                 0,                 0]];

        b = [ts(:,ptr2)+tr(:,ptr2)-s1-s2-2*K*log(lams.*lamr);...
             ts(:,ptr2)-tr(:,ptr2)-s1+s2-2*G*log(lams./lamr);...
                               sold(end,ptr1)-wold*(1./lams)];

    elseif strcmp(g_strainmeasure,'pepicelli')

        A11 = (2*K)./(lamr.*lams.^2) - (2*K*log(lamr.*lams))./(lamr.*lams.^2);
        A11 = diag(A11);
        A12 = (2*log(lamr.*lams))./(lamr.*lams);
        b1 = -(2*K*log(lamr.*lams))./(lamr.*lams);

        A21 = (2*G)./lams.^3;
        A21 = diag(A21);
        A23 = 1./lamr.^2 - 1./lams.^2;
        b2 = -G*(1./lamr.^2 - 1./lams.^2);

        A  = [[                   A11,              A12,        zeros(N,1)];...
              [                   A21,        zeros(N,1),              A23];...
              [-wold*diag(1./lams.^2),                 0,                 0]];

        b = [ts(:,ptr2)+tr(:,ptr2)-s1-s2+b1;...
             ts(:,ptr2)-tr(:,ptr2)-s1+s2+b2;...
                               sold(end,ptr1)-wold*(1./lams)];

    else

        error('strainmeasure not available!')

    end
    
    % Erase dependence and equation for lambda^r(0) were its singular
    A([1,N+1],:) = [];
    A(:,1) = [];
    b([1,N+1]) = [];
  
    u0 = A\b;
%    u0 = filterSvd(A,b,1e6);   % If ill conditioned try svd
    
    u = [0;u0];
    
% update variables    
    lams = lams + params_num.alpha_sfe*u(1:N);
    lams(1) =  -dold(1,2:end,ptr2)*lams(2:end)./dold(1,1,ptr2);
    K = K + params_num.alpha_sfe*u(N+1);
    G = G + params_num.alpha_sfe*u(N+2);
    
%   Possible to plot the evolution of K and G
    KK(iter) = K;
    GG(iter) = G;
    RS(iter) = rms(b);
    
    iter = iter+1;
    
    % reassess s in r(s) in lambda^r, sigma^s and sigma^r
    sprime = Itg*(1./lams);
    r0 = interp1(sold(:,ptr2)*sprime(end)/sold(end,ptr2),rold(:,ptr1),sprime,'spline');
    s1 = interp1(sold(:,ptr2)*sprime(end)/sold(end,ptr2),ts(:,ptr1),sprime,'spline');
    s2 = interp1(sold(:,ptr2)*sprime(end)/sold(end,ptr2),tr(:,ptr1),sprime,'spline');
    
% re-evaluate lambda^r
    lamr = rold(:,ptr2)./r0;
    lamr(1) =  -dold(1,2:end,ptr2)*lamr(2:end)./dold(1,1,ptr2);
% bring lambda^s and lambda^r together    
    dlam = lamr(1)-lams(1);
    lamr = lamr-0.25*dlam;
    lams = lams+0.25*dlam;
   
% % if the error is small and re-increases stop the iteration    
%     if (oldb<1e-1)&&(oldb<rms(b))
%         break;
%     end
    oldb = rms(b);

        if verbose
            fprintf('iter %d: rms(u(N+1:N+2)) = %d\n',iter,rms(u(N+1:N+2)));
        end
end

% print result to screen
if verbose
    display(strcat('SFE converged to K=',num2str(K),' and G=',num2str(G),' Residual=',num2str(oldb)));
end

% save fit error and also the derivative of the moduli in the error, this
% is another possible way to assess the quality of the fitting
g_error = [oldb;(KK(iter-1)-KK(iter-2))/(RS(iter-1)-RS(iter-2));(GG(iter-1)-GG(iter-2))/(RS(iter-1)-RS(iter-2))];
GK = [G, K];

end

