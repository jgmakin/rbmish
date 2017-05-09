function xpctT = enforceNoisilyObservedControls(xpctT,C,H,T)
% Somewhat recherche, perhaps.  The idea is to treat some of the states as
% "controls," which are observed but noisily.  This dynamical system is
% given by the matrices in the structure dynamics.  This function enforces
% the original system by making:
%
%   (0) correct muU and SigmaU in IC cumulants
%   (1) <Uf*Zp> = <Uf>*<Zp>                     => zeros in Gamma
%   (2) <Xf*Uf> = <Xf>*<Uf>                     => Cov[Xf,Uf|Zp] = 0
%   (3) <Y*U> = ...                             => zeros in Lambda (top)
%   (4) <V*X> = ...                             => zeros in Lambda (bottom)
%   (5) <Y*V> =                                 => Cov[Y,V|Z] = 0 
%   (6) <V*Y> =                                 => [[Cov[Y|Z] = Cov[Y|X]???
%   (&) <X0*U0> = <X0>*<U0>                     => Cov[X0,U0] = 0
%
%%% 
%
%%% <Y*U> =? <Y><U>??
%%% <V*X> =? <V><X>??
%

fprintf('enforcing noisily observed controls...\n');

% Ns
[Ny,Nx] = size(C);
[Nv,Nu] = size(H);

% muU = (<U0> + (T-1)<Uf>)/T 
muU = (xpctT.mu0(Nx+(1:Nu),1) + (T-1)*xpctT.XfXp(Nx+(1:Nu),end))/T;
xpctT.mu0(Nx+(1:Nu),1) = muU;       
xpctT.XfXp(Nx+(1:Nu),end) = muU;

% UU = (<U0*U0> + (T-1)<Uf*Uf>)/T, 
UU = (xpctT.x0x0(Nx+(1:Nu),Nx+(1:Nu)) +...
    (T-1)*xpctT.XfXf(Nx+(1:Nu),Nx+(1:Nu)))/T;
xpctT.x0x0(Nx+(1:Nu),Nx+(1:Nu)) = UU; 
xpctT.XfXf(Nx+(1:Nu),Nx+(1:Nu)) = UU;

% zeros in Gamma
Uf = xpctT.XfXp(Nx+(1:Nu),end);
Zp = xpctT.XpXp(1:end-1,end);
xpctT.XfXp(Nx+(1:Nu),1:end-1) = Uf*Zp';
%%% xpctT.XfXp/xpctT.XpXp

% zeros in SigmaZ 
XfbarUbar = xpctT.XfXp(1:Nx,end)*xpctT.XfXp(Nx+(1:Nu),end)';
xpctT.XfXf(1:Nx,Nx+(1:Nu)) = XfbarUbar;
xpctT.XfXf(Nx+(1:Nu),1:Nx) = XfbarUbar';
%%% xpctT.XfXf - xpctT.XfXp/xpctT.XpXp*xpctT.XfXp'

% zeros in Lambda
C = xpctT.YX(1:Ny,1:Nx)/xpctT.XX(1:Nx,1:Nx);
xpctT.YX(1:Ny,Nx+(1:Nu)) = C*xpctT.XX(1:Nx,Nx+(1:Nu));
H = xpctT.YX(Ny+(1:Nv),Nx+(1:Nu))/xpctT.XX(Nx+(1:Nu),Nx+(1:Nu));
xpctT.YX(Ny+(1:Nv),1:Nx) = H*xpctT.XX(Nx+(1:Nu),1:Nx);
%%% xpctT.YX/xpctT.XX

% zeros in SigmaW|Z
xpctT.YY(1:Ny,Ny+(1:Nv)) = C*xpctT.YX(Ny+(1:Nv),1:Nx)';  %%% these are 
xpctT.YY(Ny+(1:Nv),1:Ny) = H*xpctT.YX(1:Ny,Nx+(1:Nu))';  %%%   the same
%%% xpctT.YY - xpctT.YX/xpctT.XX*xpctT.YX'

% zeros in Upsilon0
X0barU0bar = xpctT.mu0(1:Nx,1)*xpctT.mu0(Nx+(1:Nu),1);
xpctT.x0x0(1:Nx,Nx+(1:Nu)) = X0barU0bar;
xpctT.x0x0(Nx+(1:Nu),1:Nx) = X0barU0bar';
%%% xpctT.x0x0 - xpctT.mu0*xpctT.mu0'



end