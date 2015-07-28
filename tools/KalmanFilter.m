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
%       <-log q(y_0,...,y_t|params)>_p(y_0,...,y_t)
%
%   The parameters of the dynamical system need to be stored in the
%   structure params.  The mean state noise (muX), emission noise (muY),
%   *input* matrix (G), inputs (U), and input variance (Q) are assumed to
%   be zero if not found as fields in the params structure.  The covariance
%   matrices of the emission and state noise can be either constant (one
%   matrix) or time varying (a Ndims x Ndims x T tensor).
%
%   The equations for the full model are:
%
%       X_{t+1} ~ N(A*X_t + G*U + muX, SigmaX)
%           Y_t = N(C*X_t + muY, SigmaY)
%           V_t = N(H*U_t + muV, SigmaV)
%
%           X_0 ~ N(x0,inv(Info0))
%
%   NB: to run with fully observed inputs, you need to change this code....

%-------------------------------------------------------------------------%
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
%   params.SigmaY is the size of a covariance matrix, it gets repmat'ed
%   accordingly
% Revised: 04/30/13
%   -functionized
%   -made it work
% Created: 04/29/13
%   by JGM
%-------------------------------------------------------------------------%

% KF params
T = params.T;
A = params.A;
C = params.C;
SigmaX = params.SigmaX;
SigmaY = params.SigmaY;

% numbers
Nstates = size(A,1);
Nemits = size(Y,1);

% optional params
if numel(SigmaY) == Nemits^2, SigmaY = repmat(SigmaY,[1,1,T]); end
if numel(SigmaX) == Nstates^2, SigmaX = repmat(SigmaX,[1,1,T]); end
if isfield(params,'muX'),muX = params.muX;else muX = zeros(1,'like',Y); end
if isfield(params,'muY'),muY = params.muY;else muY = zeros(1,'like',Y); end
if isfield(params,'B'), B = params.B; else B = ones(1,'like',Y); end
if ~isempty(varargin), U = varargin{1}; else U = zeros([1,T],'like',Y); end



% malloc
XHATMU = nan(Nstates,T,'like',Y);
CVRNMU = nan(Nstates,Nstates,T,'like',Y);
XHATTU = nan(Nstates,T,'like',Y);
INFOTU = nan(Nstates,Nstates,T,'like',Y);

% initialize with the prior params
xhatTU = params.mu0;
InfoTU = params.Info0;

% cross entropy that accumulates across all t (times 2!) 
XNtrpY = T*Nstates*log(2*pi);
for t = 1:T

    % store (xhat_{t|t-1},P_{t|t-1})
    XHATTU(:,t) = xhatTU;
    INFOTU(:,:,t) = InfoTU; 

    % measurement update
    if isinf(det(SigmaY(:,:,t))),InfoYX = 0; else InfoYX = inv(SigmaY(:,:,t)); end
    %%%%
%     if rank(C'*InfoYX*C + InfoTU) < size(C'*InfoYX*C + InfoTU,1)
%         CvrnMU = params.Info0;
%         CvrnMU(1) = 0.1218;
%         fprintf('doing some bad things at time %i...\n',t);
%     else
%         CvrnMU = inv(C'*InfoYX*C + InfoTU);
%     end
    %%%%
    % CvrnMU = [0.5401 -0.8408; -0.8408 1.3096];
    CvrnMU = inv(C'*InfoYX*C + InfoTU);
    xhatMU = CvrnMU*(C'*InfoYX*(Y(:,t)-muY) + InfoTU*xhatTU);   
    
    % store (xhat_{t|t},P_{t|t})
    XHATMU(:,t) = xhatMU;
    CVRNMU(:,:,t) = CvrnMU;
    
    % "cross-entropy"
    innov = Y(:,t) - muY - C*xhatTU;
    InfoY = InfoYX - InfoYX*C*CvrnMU*C'*InfoYX;
    XNtrpY = XNtrpY - log(det(InfoY)) + innov'*InfoY*innov;
   
    % time update
    InfoTU = inv(A*CvrnMU*A' + SigmaX(:,:,t));
    xhatTU = A*xhatMU + B*U(:,t) + muX;

end
XNtrpY = XNtrpY/2;      % common factor of 1/2 left out above
XNtrpY = XNtrpY/T;      % i.e., time-averaged cross entropy (sensible)


% collect outputs
KFdstrbs.XHATMU = XHATMU;
KFdstrbs.CVRNMU = CVRNMU;
KFdstrbs.XHATTU = XHATTU;
KFdstrbs.INFOTU = INFOTU;
KFdstrbs.XNtrpY = XNtrpY;


% announce
%%% fprintf('\n\nFiltered...\n\n');


end




