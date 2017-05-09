function [B,Rsq,Yres,XpctTY] = IRLS(X,Y,dstrb)
% IRLS  Iteratively reweighted least squares
%
% USAGE:
%   B = IRLS(X,Y,dsrtb);
%
% For each dstrb, IRLS.m assumes the canonical response function:
%   'Bernoulli'     logistic
%   'Poisson'       exponential
%   ...
%
% Both X and Y can be matrices, whose *rows* are "samples" or "trials."
% Hence:
%
%   size(X)     = Nsamples x Nfeatues
%   size(Y)     = Nsamples x Noutputvars
%   size(beta)  = Nfeatures x Noutputvars
%
% NB: Does *not* include the column of ones by dafault.

%-------------------------------------------------------------------------%
% Revised; 02/13/17
%   -added a little recursive function to deal with data too big to
%   compute all outer products of at once.
% Created: 02/11/17
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Nsamples,Nx] = size(X);
Ny = size(Y,2);
NiterMax = 100;

% init
thr = 0.0005;
B = rand(Nx,Ny,'like',Y)/1000;
XpctTY = zeros(size(Y),'like',Y);

% sufficient stats that can be calculate once for all
XY = X'*Y;                                      % (Nx x Ny)
if Nsamples*Nx^2 < 5e8
    XX = permute(X,[2,3,1]).*permute(X,[3,2,1]);
else
    XX = [];
end
%%% but NB: *these* are *not* summed!! but rather (Nx x Nx x Nsamples)

fprintf('IRLS loop...\n')
% the outputs are independent, but unfortunately this cannot be vectorized
for iy = 1:Ny
    
    Niter = 0;
    db = NaN(Nx,1,'like',Y);
    b = B(:,iy);
    
    while ~all(abs(gather(db./b))<thr)&&(Niter<NiterMax)
        
        % get gradient and (negative) Hessian
        [negH,grad,XpctTy] = getGradientAndHessian(X,XX,XY(:,iy),b,dstrb);
        
        % Newton-Raphson: deltaphi = -inv(d^2L/dphi^2)*(dL/dphi)
        db = negH\grad;        
        b = b + db;
        Niter = Niter+1;
        fprintf('.');
    end
    
    B(:,iy) = b;
    XpctTY(:,iy) = XpctTy;
    fprintf('\n');
    
end

% coefficient of determination
Yres = Y - XpctTY;
SStot = sum((Y - mean(Y,1)).^2);
SSerr = sum(Yres.^2,1);
Rsq = 1 - SSerr./SStot;


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [negH,grad,XpctTY] = getGradientAndHessian(X,XX,Xy,beta,dstrb)

% E[T(Y)]   (Nsamples x 1)
XpctTY = invParamMap(X,beta,zeros(size(X,1),1,'like',X),{dstrb},[],[]);

% Var[T(Y)] (Nsamples x 1)
switch dstrb
    case 'Bernoulli'
        VrncTY = XpctTY.*(1 - XpctTY);
    case 'StandardNormal'
        VrncTY = ones(size(XpctTY),'like',XpctTY);
    case 'Poisson'
        
    case 'Gamma'
      
    otherwise
        error('unrecognized distribution for IRLS -- jgm');
        
end


% gradient and (negative) Hessian
grad = (Xy - X'*XpctTY);                                    % (Nx x 1)
if isempty(XX)
    % in case it was impossible to store all the outer products
    negH = weightedSumOuterProducts(X,VrncTY);
else
    negH = sum(XX.*permute(VrncTY,[2,3,1]),3);% (Nx x Nx)
end






end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function XwX = weightedSumOuterProducts(X,w)
% The variables are too big, but don't loop (the other extreme) through
% outer products: break the data in half recursively until the pieces are
% small enough to fit in memory

[Nsamples,Nx] = size(X);

if Nsamples*Nx^2 < 5e7
    XX = permute(X,[2,3,1]).*permute(X,[3,2,1]);
    XwX = sum(XX.*permute(w,[2,3,1]),3);
else
    midind = floor(size(X,1)/2);
    XwX = weightedSumOuterProducts(X(1:midind,:),w(1:midind));
    XwX = XwX + weightedSumOuterProducts(X((midind+1):end,:),w((midind+1):end));
    
end

end
%-------------------------------------------------------------------------%




























