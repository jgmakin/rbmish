function RTSSdstrbs = RTSsmoother(params,KFdstrbs)
% RTSsmoother   Rauch-Tung-Striebel smoother
%
% RTSSdstrbs = RTSsmoother(LDSparams,KFdstrbs) returns the cumulants of 
% the (Gaussian) posterior distribution over the hidden state:
%
%   XHAT     = E[X_t|y_0,...,y_T]
%   XCVRN    = Cov[X_t|y_0,...,y_T]
%   XfXpCVRN = Cov[X_{t+1},X_t|y_0,...,y_T]
%
% for all t.  These three variables are returned as fields in the output 
% structure, RTSSdstrbs.  The input structure LDSparams must contain the
% state transition matrix A as a field.  KFdstrbs is the structure returned
% by KalmanFilter.m.

%-------------------------------------------------------------------------%
% Revised: 09/22/16
%   -added a version (almost entirely redundant) to work with CVRNTU rather
%   than INFOTU, in accordance with changes made to KalmanFilter.m
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


if isfield(KFdstrbs,'CVRNTU')
    [XHAT,XCVRN,XfXpCVRN] = RTSSv1(params.A,KFdstrbs.XHATMU,...
        KFdstrbs.CVRNMU,KFdstrbs.XHATTU,KFdstrbs.CVRNTU);
else
    [XHAT,XCVRN,XfXpCVRN] = RTSSv2(params.A,KFdstrbs.XHATMU,...
        KFdstrbs.CVRNMU,KFdstrbs.XHATTU,KFdstrbs.INFOTU);
    
end

% announce
%%% fprintf('\n\nSmoothed...\n\n');

% collect outputs
RTSSdstrbs.XHAT = XHAT;
RTSSdstrbs.XCVRN = XCVRN;
RTSSdstrbs.XfXpCVRN = XfXpCVRN;

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [XHAT,XCVRN,XfXpCVRN] = RTSSv1(A,XHAT,XCVRN,XHATTU,XfXpCVRN)
% NB that XHATMU, CVRNMU, and CVRNTU are just initialized as XHAT, XCVRN,
% and XfXpCVRN.  The first two are only right for the final (T) sample, and
% the last not right at all!  But the loop below updates all other samples.
% This was done to save memory.


% Ns
T = size(XHAT,2);

% initialize
Ptplus1t = XfXpCVRN(:,:,T);                     % i.e, CVRNTU(:,:,T)
XfXpCVRN = XfXpCVRN(:,:,1:T-1);                 % only need T-1 frames

% backwards recursion (see labnotes)
for t = (T-1):-1:1
    
    % useful quantities
    INNOV = XHAT(:,t+1) - XHATTU(:,t+1);
    Pt = XCVRN(:,:,t);                          % pre-updated XCVRN=CVRNMU
    Jt = Pt*A'/Ptplus1t;
    Ptplus1t = XfXpCVRN(:,:,t);                 % i.e, CVRNTU(:,:,t)
    
    % Cov[X_t,X_{t+1}|y], E[X_t|y], Cov[X_t|y]
    XfXpCVRN(:,:,t) = XCVRN(:,:,t+1)*Jt';
    XHAT(:,t) = XHAT(:,t) + Jt*INNOV;           % pre-updated XHAT=XHATMU
    XCVRN(:,:,t) = Pt + Jt*(XfXpCVRN(:,:,t) - A*Pt);
        
    %%%% equivalent, but ever so slightly slower %%%%
    % cvrnInnov = XCVRN(:,:,t+1) - CVRNTU(:,:,t+1);
    % XCVRN(:,:,t) = Pt + Jt*cvrnInnov*Jt';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [XHAT,XCVRN,XfXpCVRN] = RTSSv2(A,XHAT,XCVRN,XHATTU,XfXpCVRN)
% NB that XHATMU, CVRNMU, and INFOTU are just initialized as XHAT, XCVRN,
% and XfXpCVRN.  The first two are only right for the final (T) sample, and
% the last not right at all!  But the loop below updates all other samples.
% This was done to save memory.


% Ns
T = size(XHAT,2);

% initialize
Ptplus1tInv = XfXpCVRN(:,:,T);                  % i.e., INFOTU(:,:,T)
XfXpCVRN = XfXpCVRN(:,:,1:T-1);                 % only need T-1 frames

% backwards recursion (see labnotes)
for t = (T-1):-1:1
    
    % useful quantities
    INNOV = XHAT(:,t+1) - XHATTU(:,t+1);
    Pt = XCVRN(:,:,t);                          % pre-updated XCVRN=CVRNMU
    Jt = Pt*A'*Ptplus1tInv;
    Ptplus1tInv = XfXpCVRN(:,:,t);              % i.e, INFOTU(:,:,t)
    
    % Cov[X_t,X_{t+1}|y], E[X_t|y], Cov[X_t|y]
    XfXpCVRN(:,:,t) = XCVRN(:,:,t+1)*Jt';
    XHAT(:,t) = XHAT(:,t) + Jt*INNOV;           % pre-updated XHAT=XHATMU
    XCVRN(:,:,t) = Pt + Jt*(XfXpCVRN(:,:,t) - A*Pt);
    
end

end
%-------------------------------------------------------------------------%
