function LDSparams = learnfullyobservedLDS(LDSdata,varargin)
% Just like it says: This version of the parameter estimator learns them
% based on a fully observed model; i.e., it has access to the full
% underlying state X.
%
% Right now it learns using only one data set, i.e., 40 trajectories of
% 1000 samples apiece.


%%%%%%%%%%%%%%
% NB: you need a flag to cover the following cases:
%   (1) SigmaY is time varying (!?)
%   (2) regression with LOOCV
%   (3) regression with regularization (but that's in linrgsLOOCV, too)
%   (4) unbiased estimators?? (N-1 vs. N)
% These will be a pain in the ass, b/c getMLparams is called in EM4KF.m as
% well, and there it's a little more confusing how to compute
% cross-validated parameters (you're working with expected sufficient stats
% rather than actually samples of X...)
%%%%%%%%%%%%%%


%-------------------------------------------------------------------------%
% Revised: 08/26/14
%   -eliminated "noisily observed control"
% Revised: 01/15/14
%   -rearranged, rationalized, made more similar to EM4KF
%   -added t.muU0
% Revised: 01/06/14
%   -changed input to structure to accommodate controls
%   -added if statements for controlled LDS (necessitated by unfortunate
%   facts about cat'ing empty tensors with nonzero sizes...)
% Revised: 12/18/13
%   -put X into "canonical structure" (a dimension for Nmods)
% Revised: 11/26/13
%   -now assumes non-zero-mean noise terms (the "columns of ones") for both
%   the emissions and the transitions
% Cribbed: 11/21/13
%   from yr getLDSparams
%   by JGM
%-------------------------------------------------------------------------%


% unload structure
Z = LDSdata.Z;
Y = LDSdata.Y;

% params
[Ncases, Nstates, T] = size(Z);
Nsamples = Ncases*T;
NN = Nsamples - Ncases;
LDSparams.T = T;


% sift out data for regressions (Z is a tensor: Ncases x Nstates x T)
x00 = Z(:,:,1)';
Xp = [longdata(Z(:,:,1:end-1))'; ones(1,NN)];
Xf = longdata(Z(:,:,2:end))';
X = [longdata(Z)'; ones(1,Nsamples)];
Ybroad = longdata(Y)';

%%%%%
% adjust for (visible) control
%%%%%

% COMPUTE SUFFICIENT STATISTICS
% first-order
t.mu0 = mean(x00,2);
%%% this will require a lot more trials (than 40) to get right

% second-order: expected outer products of...
t.x0x0 = x00*x00'/Ncases;                       % ...the first
t.XpXp = Xp*Xp'/NN;                             % ...all but the last
t.XfXp = Xf*Xp'/NN;                             % ...consecutive states
t.XfXf = Xf*Xf'/NN;                             % ...all but the first
t.XX = X*X'/Nsamples;                           % ...all states
t.YX = Ybroad*X'/Nsamples;                      % ...states and emissions
t.YY = Ybroad*Ybroad'/Nsamples;                 % ...all emissions

% tell getMLparams to expect non-zero-mean noise terms (since everything is
% observed, these can't be redundant with each other...)
LDSparams.muY = [];      LDSparams.muX = [];

%%%%%
if ~isempty(varargin)
    Nx = varargin{1}; Ny = varargin{2}; Nu = varargin{3}; Nv = varargin{4};
    dynamics.C = zeros(Ny,Nx);  dynamics.H = zeros(Nv,Nu);  dynamics.T = T;
    t = enforceNoisilyObservedControls(t,dynamics);
end
%%%%%

% solve for the maximum-likelihood parameters (may be biased!)
LDSparams = getMLparams(t,LDSparams);



end

