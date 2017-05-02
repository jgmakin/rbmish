function KFdstrbs = KalmanFilter(params,Y,varargin)
% KalmanFilter  The Kalman filter
%
% USAGE:
%   KFdstrbs = KalmanFilter(params,Y);
%   KFdstrbs = KalmanFilter(params,Y,U);
%
%   The Kalman filter is the recursive solution to an inference problem;
%   sc., computing the posterior distribution:
%
%       p(x_t|y_0,...,y_t)
%
%   over the hidden state x_t, for all t, in a system with linear-Gaussian
%   emissions and linear-Gaussian dynamics.  If the initial state is
%   normally distributed, the posterior will be as well, for all t, so that
%   the filter needs only to update two cumulants, the mean (XHATMU) and
%   covariance matrix (CVRNMU).
%
%   When learning the parameters of the filter, it is also useful to know
%   the (negative) average log-likelihood of the model parameters under the
%   data, i.e. the cross entropy (XNtrp) of the model under the data:
%
%       <-log q(y_0,...,y_t|params)>_p(y_0,...,y_t).
%
%   The parameters of the dynamical system need to be stored in the
%   structure params.  The mean state noise (muX), emission noise (muYX),
%   input matrix (B), inputs (U), and input variance (Q) are assumed to
%   be zero if not found as fields in the params structure.  The covariance
%   matrices of the emission and state noise can be either constant (one
%   matrix) or time varying (a Ndims x Ndims x T tensor).
%
%   The equations for the full model are:
%
%           X_0 ~ N(x0,inv(Info0))            
%       X_{t+1} ~ N(A*X_t + B*U + muX, SigmaX)
%           Y_t ~ N(C*X_t + muYX, SigmaYX)


%-------------------------------------------------------------------------%
% Revised: 09/21/16
%   -implemented the other version of the KF (with a Kalman gain), because
%   this should be more efficient when there are more states than
%   observations---and because it's easier to accommodate zero-information
%   observations (as in the 'LTI-PPC' when the visibles are Bernoulli.)
% Revised: 04/11/16
%   -allowed for InfoYX to be passed directly as a field of params, rather
%   than expecting only SigmaYX.
%   -commented out code to invert SigmaYX when it has less than full rank or
%   infinite elements.  This may have been a mistake.
% Revised: 08/28/14
%   -changed to assume that controls are "fully observed."
%   -eliminated "adjusters"
%   -changed varargin to hold the controls, U
% Revised: 01/14/14
%   -changed outputs to single structure
%   -added output fields for UHAT, CVRNUHAT, XHATTU, and CVRNTU
% Revised: 01/09/14
%   -fixed the assignment of optional parameters to work both with and
%   without an input
% Revised: 01/06/14
%   -added funtionality for a controlled LDS, where the controls are
%   (possibly) observed noisily.
% Revised: 07/02/13
%   -added varargin for a function to "adjust" xhatTU w/i the filter loop
%   -forced filter loop to use a different emission cvrn on every trial; if
%   params.SigmaYX is the size of a covariance matrix, it gets repmat'ed
%   accordingly
% Revised: 04/30/13
%   -functionized
%   -made it work
% Created: 04/29/13
%   by JGM
%-------------------------------------------------------------------------%

% KF params
A = params.A;
C = params.C;
SigmaX = params.SigmaX;

% Ns
Nstates = size(A,1);
[Nobsvs,T] = size(Y);

% optional params
if isfield(params,'muX'),muX = params.muX;else muX = zeros(1,'like',Y); end
if isfield(params,'muYX'),muYX = params.muYX;else muYX = zeros(1,'like',Y); end
if isfield(params,'B'), B = params.B; else B = ones(1,'like',Y); end
U = defaulter('controls',zeros([1,T],'like',Y),varargin{:});
LIGHTWEIGHT = defaulter('lightweight',0,varargin{:}); % don't save all

% might as well do this all at once
Y   = Y - muYX;
muX = B*U + muX;

% Kalman filtering (two different ways)
if (Nstates > Nobsvs)&&~any(diag(params.Info0)==0)
    
    % ...then it's better to invert emission covariances
    if isfield(params,'InfoYX')
        CvrnYX  = arrayfun(@(ii)(yrinv(params.InfoYX(:,:,ii))),...
            1:size(params.InfoYX,3),'UniformOutput',false);
        CvrnYX  = cat(3,CvrnYX{:});
    else
        CvrnYX  = params.SigmaYX;
    end
    KFdstrbs = KFv1(Y,C,CvrnYX,A,muX,SigmaX,params.mu0,...
        yrinv(params.Info0),LIGHTWEIGHT);
    
else
    
    % ...it's better to invert state covariances
    if isfield(params,'InfoYX')
        InfoYX  = params.InfoYX;
    else
        InfoYX = arrayfun(@(ii)(yrinv(params.SigmaYX(:,:,ii))),...
            1:size(params.SigmaYX,3),'UniformOutput',false);
        InfoYX  = cat(3,InfoYX{:});
    end
    CtrInfo = tensorOp(repmat(C',[1,1,size(InfoYX,3)]),InfoYX);
    KFdstrbs = KFv2(Y,C,InfoYX,CtrInfo,A,muX,SigmaX,params.mu0,...
        params.Info0,LIGHTWEIGHT);
end

% announce
%%% fprintf('\n\nFiltered...\n\n');


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function KFdstrbs = KFv1(Y,C,CvrnYX,A,muX,SigmaX,xhatTU,CvrnTU,LIGHTWEIGHT)
% This is the standard version you see in textbooks

% Ns
Nstates = size(SigmaX,1);
T = size(Y,2);

% malloc
XNtrpY = 0;
XHATMU = nan(Nstates,T,'like',Y);
if ~LIGHTWEIGHT
    CVRNMU = nan(Nstates,Nstates,T,'like',Y);
    XHATTU = nan(Nstates,T,'like',Y);
    CVRNTU = nan(Nstates,Nstates,T,'like',Y);
end

Nsteps = T;
for t = 1:T
    
    % store (xhat_{t|t-1},P_{t|t-1})
    if ~LIGHTWEIGHT, XHATTU(:,t) = xhatTU; CVRNTU(:,:,t) = CvrnTU; end
    
    % measurement update and cross-entropy calculation
    thisCvrnYX = CvrnYX(:,:,min(t,size(CvrnYX,3)));
    if isinf(sum(sum(thisCvrnYX))) % no information at this step!
        %%% seems wrong--what if there's info about **some** obsvs?
        
        % measurement update
        xhatMU = xhatTU;
        CvrnMU = CvrnTU;
        
        % no cross-entropy contribution--decrement the denominator
        Nsteps = Nsteps-1;
    else
        
        % measurement update
        K = CvrnTU*C'/(thisCvrnYX + C*CvrnTU*C');
        innov = Y(:,t) - C*xhatTU;
        xhatMU = xhatTU + K*innov;
        CvrnMU = CvrnTU - K*C*CvrnTU;
        
        % cross-entropy of Y under the model: <-log{q(y;theta)}>_p(y)
        CvrnY = C*CvrnTU*C' + thisCvrnYX;
        thisXNtrpY = GaussianCrossEntropy(innov,CvrnY,'cvrn',0);
        XNtrpY = XNtrpY + thisXNtrpY;
    end
    
    % store (xhat_{t|t},P_{t|t})
    XHATMU(:,t) = xhatMU;
    if ~LIGHTWEIGHT, CVRNMU(:,:,t) = CvrnMU; end
    
    % time update
    CvrnTU = A*CvrnMU*A' + SigmaX(:,:,min(t,size(SigmaX,3)));
    xhatTU = A*xhatMU + muX(:,t);
    
end
XNtrpY = XNtrpY/Nsteps;     % i.e., time-averaged cross entropy (sensible)

% collect outputs
KFdstrbs.XHATMU = XHATMU;
if ~LIGHTWEIGHT
    KFdstrbs.CVRNMU = CVRNMU;
    KFdstrbs.XHATTU = XHATTU;
    KFdstrbs.CVRNTU = CVRNTU;
end
KFdstrbs.XNtrpY = XNtrpY;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function KFdstrbs = KFv2(Y,C,InfoYX,CtrInfo,A,muX,SigmaX,xhatTU,InfoTU,LIGHTWEIGHT)
% This is the version that looks more like "multisensory integration"

% Ns
Nstates = size(SigmaX,1);
T = size(Y,2);

% malloc
XNtrpY = 0;
XHATMU = nan(Nstates,T,'like',Y);
% do you really need these? They can take up a lot of space
if ~LIGHTWEIGHT
    CVRNMU = nan(Nstates,Nstates,T,'like',Y);
    XHATTU = nan(Nstates,T,'like',Y);
    INFOTU = nan(Nstates,Nstates,T,'like',Y);
end

% step through time
for t = 1:T
    
    % store (xhat_{t|t-1},P_{t|t-1})
    if ~LIGHTWEIGHT, XHATTU(:,t) = xhatTU; INFOTU(:,:,t) = InfoTU; end
    
    % measurement update
    thisCtrInfo = CtrInfo(:,:,min(t,size(CtrInfo,3)));
    InfoMU = thisCtrInfo*C + InfoTU;
    xhatMU = InfoMU\(thisCtrInfo*Y(:,t) + InfoTU*xhatTU);
    
    % store (xhat_{t|t},P_{t|t})
    XHATMU(:,t) = xhatMU;
    if ~LIGHTWEIGHT, CVRNMU(:,:,t) = inv(InfoMU); end
    
    % cross-entropy of Y under the model: <-log{q(y;theta)}>_p(y)
    innov = Y(:,t) - C*xhatTU;
    thisInfoYX = InfoYX(:,:,min(t,size(InfoYX,3)));
    InfoY = thisInfoYX - thisCtrInfo'/InfoMU*thisCtrInfo;
    XNtrpY = XNtrpY + GaussianCrossEntropy(innov,InfoY,'info',0);
    
    % time update
    InfoTU = inv(A/InfoMU*A' + SigmaX(:,:,min(t,size(SigmaX,3))));
    xhatTU = A*xhatMU + muX(:,t);
    
end
XNtrpY = XNtrpY/T;      % i.e., time-averaged cross entropy (sensible)

% collect outputs
KFdstrbs.XHATMU = XHATMU;
if ~LIGHTWEIGHT
    KFdstrbs.CVRNMU = CVRNMU;
    KFdstrbs.XHATTU = XHATTU;
    KFdstrbs.INFOTU = INFOTU;
end
KFdstrbs.XNtrpY = XNtrpY;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function Minv = yrinv(M)
% This is somewhat dubious.  The idea is that you often have diagonal
% covariance matrices with an Inf or a zero on the diagonal, representing
% absolute certainty or uncertainty.  Using pinv to invert this replaces
% zeros with zeros, which would turn (e.g.) absolute uncertainty (zero
% information) into absolute certainty (zero covariance)---a mistake.


switch det(M)
    case {inf,0}
        if isdiag(M)
            Minv = diag(1./diag(M));
        else
            Minv = pinv(M);
        end
    otherwise
        Minv = inv(M);
end

end
%-------------------------------------------------------------------------%
