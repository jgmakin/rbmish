function RTSSdstrbs = RTSsmoother(params,KFdstrbs)
% Returns the cumulants of the (Gaussian) posterior distribution p(x_y|y), 
% where y is *all* the observations.
%
% Also returns the posterior covariance Cov[x_t,x_{t+1}|y] =: S.

%-------------------------------------------------------------------------%
% Revised: 08/28/14
%   -removed update equations for "noisily observed control," since this
%   can just be treated as a state (with holes in the parameter matrices)
% Revised: 01/20/14
%   -added updates for the control
% Revised: 01/14/14
%   -rewrote in terms of the other set of update equations (see lab notes)
% Revised: 01/13/14
%   -added accommodations for (noisily observed) inputs to the LDS
% Revised: 10/28/13
%   -compacted initialization and malloc
% Created: 10/23/13
%   by JGM
%-------------------------------------------------------------------------%

% filtered distributions
XHATMU = KFdstrbs.XHATMU;
CVRNMU = KFdstrbs.CVRNMU;
XHATTU = KFdstrbs.XHATTU;
INFOTU = KFdstrbs.INFOTU;

% params
T = size(XHATMU,2);
A = params.A;
Nstates = size(A,1);

% initialize with the final filter distribution/malloc
XHAT = XHATMU;
XCVRN = CVRNMU;
XfXpCVRN = zeros(Nstates,Nstates,T-1,'like',XHAT);


% backwards recursion (see labnotes)    
for t = (T-1):-1:1
    
    % useful quantities    
    INNOV = XHAT(:,t+1) - XHATTU(:,t+1);
    Pt = CVRNMU(:,:,t);
    Jt = Pt*A'*INFOTU(:,:,t+1);

    % Cov[X_t,X_{t+1}|y], E[X_t|y], Cov[X_t|y]
    XfXpCVRN(:,:,t) = XCVRN(:,:,t+1)*Jt';
    XHAT(:,t) = XHATMU(:,t) + Jt*INNOV;
    XCVRN(:,:,t) = Pt + Jt*(XfXpCVRN(:,:,t) - A*Pt);
    %%%% equivalent, but slower %%%%
    % cvrnInnov = XCVRN(:,:,t+1) - inv(INFOTU(:,:,t+1));
    % XCVRN(:,:,t) = Pt + Jt*cvrnInnov*Jt';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% announce
%%% fprintf('\n\nSmoothed...\n\n');

% collect outputs
RTSSdstrbs.XHAT = XHAT;
RTSSdstrbs.XCVRN = XCVRN;
RTSSdstrbs.XfXpCVRN = XfXpCVRN;

end


