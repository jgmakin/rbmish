function LDSparams = learnfullyobservedLDS(Y,X,varargin)
% Just like it says: This version of the parameter estimator learns them
% based on a fully observed model; i.e., it has access to the full
% underlying state X.
% 
% USAGE:
%    LDSparams = learnfullyobservedLDS(Y,X);
%


%%%%%%%%%%%%%%
% NB: you need a flag to cover the following cases:
%   (1) SigmaYX is time varying (!?)
%   (2) regression with LOOCV
%   (3) regression with regularization (but that's in linregress.m, too)
%   (4) unbiased estimators?? (N-1 vs. N)
% These will be a pain in the ass, b/c ML4LDS is called in EM4KF.m as
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

% params
[Ncases, Nstates, T] = size(X);
Nsamples = Ncases*T;
NN = Nsamples - Ncases;
LDSparams.T = T;
altdims = defaulter('alternative dimensions',[],varargin{:});


% sift out data for regressions (X is a tensor: Ncases x Nstates x T)
x00 = X(:,:,1)';
Xp = [longdata(X(:,:,1:end-1))'; ones(1,NN)];
Xf = longdata(X(:,:,2:end))';
Xbroad = [longdata(X)'; ones(1,Nsamples)];
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
t.XX = Xbroad*Xbroad'/Nsamples;                 % ...all states
t.YX = Ybroad*Xbroad'/Nsamples;                 % ...states and emissions
t.YY = Ybroad*Ybroad'/Nsamples;                 % ...all emissions

% tell ML4LDS to expect non-zero-mean noise terms (since everything is
% observed, these can't be redundant with each other...)
LDSparams.muYX = [];            LDSparams.muX = [];

if ~isempty(altdims)
    Nx = altdims(1); Ny = altdims(2); Nu = altdims(3); Nv = altdims(4);
    t = enforceNoisilyObservedControls(t,zeros(Ny,Nx),zeros(Nv,Nu),T);
end


% solve for the maximum-likelihood parameters (may be biased!)
LDSparams = ML4LDS(t,LDSparams,varargin{:});

end