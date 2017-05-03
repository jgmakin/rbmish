function [X,Q] = getLatentsLTI(Nexamples,T,yrclass,source,params,varargin)
% getLatentsLTI
%
% USAGE:
%   [X,Q] = getLatentsLTI(Nexamples,T,yrclass,params.NS,params);
%
%   [X,Q] = getLatentsLTI(Nexamples,T,yrclass,params.NS,'gains',G);
%
% Generate "trajectories" through latent space.  The tensor X has
% size (Nexamples x Ndims x Nltnts), where Nltnts is typically just 1: it 
% corresponds to the state space (which will itself have Ndims > 1).  To 
% convert from X to *stimuli*--e.g., for encoding in PPCs, or for computing
% decoding errors--use the conversion functions stored in the cell array 
% Q.latent2stim.


%-------------------------------------------------------------------------%
% Revised: 01/05/17
%   -changed output from "stimuli" S to latent states X.  Formerly, both
%   were always carried around (the latter as a field in Q)--but S can
%   always be recovered from X!--given the appropriate conversion function.
%   Those are now stored in the cell array Q.latent2stim.
% Revised: 08/26/16 (JGM)
%   -rationalized the use of smin and smax: now they refer exclusively to
%   the stimulus S; to get the limits of the *state* Z, one uses the pinv
%   of params.dynamics.C.  NB that this will usually give meaningless
%   answers for the velocity limits---but you don't use those.
% Revised: 08/25/16 (JGM)
%   -combined cases (rEFH for 1DrEFH, 1DrEFHwithEC, 2DrEFH, 2DrEFHwithEC)
%   -moved noise generation out of loop
% Revised: 09/03/14 (JGM)
%   -added cases for wrapping + control
% Revised: 12/18/13 (JGM)
%   -added case "controlled"
% Revised: 11/20/13 (BKD)
%   -changed "resetting" condition to use wall_reset.m
%   -fixed bug that filterdata.states wasn't getting updated after resets
% Revised: 07/01/13 (JGM)
%   -renamed "flag" to "RESTART," since "flag" is protected
%   -changed "RESTART" from a huge tensor to a list of restart indices
%   -put the "new" outputs, "RESTART" and and "states," into a structure,
%   "filterdata"
% Revised: 06/27/13 (BKD)
%   -WARNING: OUTPUTS WERE CHANGED HERE
%   -added 'restart' condition
%   -now outputs states, which contains velocity
% Revised: 05/11/13
%   -wrote a function to generate ICs from either a normal or a uniform
%   prior, according to your conventions
% Created: 05/07/13
%   by JGM
%-------------------------------------------------------------------------%

% TO DO:
% (1) allow stims directly as varargin, as in getStimuliMultisensory...

% params
Ndims   = params.Ndims;
mods    = params.mods;
smin    = params.smin;
smax    = params.smax;
gmin    = params.gmin;
gmax    = params.gmax;
dynamics= params.dynamics;

% init
TOPLOT = 0;
T = defaulter('sequencelength',T,varargin{:});
pG0.mu  = NaN(length(gmin),1,yrclass);
pG0.cov = Inf;              % => uniform distribution
Ntraj   = floor(Nexamples/T);

% unload dynamical parameters
A       = dynamics.A;
C       = dynamics.C;
SigmaX  = dynamics.SigmaX;
muX0    = dynamics.muX0;
SigmaX0 = dynamics.SigmaX0;
muV0    = dynamics.muV0;
SigmaV0 = dynamics.SigmaV0;
if isfield(dynamics,'muX'),
    muX = dynamics.muX;
else
    muX = zeros(size(A,1),1);
end
zmin = pinv(C)*smin(:);     % NB: for noninvertible M, this
zmax = pinv(C)*smax(:);     %   can produce bad results!!

% use a GPU?
if strcmp(yrclass,'gpuArray')
    A   = gpuArray(A);
    muX = gpuArray(muX);            SigmaX  = gpuArray(SigmaX);
    muX0 = gpuArray(muX0);          SigmaX0 = gpuArray(SigmaX0);
    muV0 = gpuArray(muV0);          SigmaV0 = gpuArray(SigmaV0);
    zmin = gpuArray(zmin);          zmax    = gpuArray(zmax);
end


% initial conditions
mrgn = 0.05; %%% really shouldn't be expressed in absolute units...
x0 = UniformNormalDiracSampler(muX0,SigmaX0,Ntraj,zmin(1:Ndims),zmax(1:Ndims),mrgn);
v0 = UniformNormalDiracSampler(muV0,SigmaV0,Ntraj,[],[],mrgn);
if any(strcmp(mods,'Efference-Copy'))
    muU0    = dynamics.muU0;
    SigmaU0 = dynamics.SigmaU0;
    umin    = smin(:,strcmp(mods,'Efference-Copy'));
    umax    = smax(:,strcmp(mods,'Efference-Copy'));
    u0 = UniformNormalDiracSampler(muU0,SigmaU0,Ntraj,umin,umax,mrgn);
else
    u0 = [];
end

% init
X = zeros(Ntraj,size(A,1),T,yrclass);
X(:,:,1) = [x0 v0 u0];
W = shortdata(Ntraj,3,mvnrnd(zeros(Ntraj*T,length(muX),yrclass),SigmaX) + muX');
 
% the real work of this function: simulate the LTI system:
% X_{t+1} = A*X_t + W_t
for t = 1:(T-1), X(:,:,t+1) = X(:,:,t)*A' + W(:,:,t); end

% for debugging purposes
if TOPLOT&&(Ndims==2)
    NSind = strcmp(mods,source);
    plotPPCtraj(X,S,smin(:,NSind),smax(:,NSind));
end

% paralleling static data-generation data structures
X = longdata(X);

% build the functions that convert from latent vars to "stimuli"
Q.latent2stim = setConversionFxns(mods,source,C,params);

% draw the gains
Q.G = defaulter('gains',UniformNormalDiracSampler(pG0.mu,pG0.cov,...
    Ntraj*T,gmin,gmax,0),varargin{:});

% other useful things to store
Q.restarts  = 1:T:(Ntraj*T);
Q.T         = T;

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function latent2stim = setConversionFxns(mods,source,C,params)
% Note that X officially has size (Nexamples x Ndims x Nltnts) (cf.
% getLatentsMultisensoryIntegration.m), but that Nltnts = 1: the state
% space counts as a single latent variable with dimension Ndims.

% init
Nmods = length(mods);
Ndims = params.Ndims;
if isfield(params,'roboparams')
    roboparams = params.roboparams;
else
    roboparams = [];
end

% malloc
latent2stim = cell(1,Nmods);

% now store the converting functions, for use *on demand*
for iMod = 1:Nmods
    
    % get the appropriate transformation function
    outmod = mods{iMod};
    f = setSensoryTransformations(outmod,source,roboparams);
    
    % now make sure the right data go into it (complicated)
    outselector = vect(repmat(strcmp(mods,outmod),[Ndims,1]));
    M = outselector'*C;
    latent2stim{iMod} = @(X)(f(X*M',0));
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotPPCtraj(Z,src,thmin,thmax)

th0 = get2DOutline(thmin,thmax,42);
for iCase = 1:size(Z,1)
    
    % plot
    figure(1673); clf; hold on;
    
    plot(squeeze(Z(iCase,1,:)),squeeze(Z(iCase,2,:)),'r');
    plot(squeeze(src(iCase,1,1,:)),squeeze(src(iCase,2,1,:)),'m');
    plot(th0(:,1),th0(:,2),'k');
    
    axis equal;
    title(num2str(iCase));
    hold off;
    pause()
end

end
%-------------------------------------------------------------------------%