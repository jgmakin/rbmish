function [S,Q] = getLatentsHHS(Nexamples,T,yrclass,source,params,varargin)


%-------------------------------------------------------------------------%
% Revised: 12/31/16
% Cribbed: 12/31/16
%   -from trajectoryGen.m
%   by JGM
%-------------------------------------------------------------------------%

% TO DO:
% (1) find the files to load ('HHS/extracteddata/KFtuningdataHHSD080620')
% (2) is "source" ever used??
% (3) Is there any point in keeping C and H separate?

% params
Ndims   = params.Ndims;
mods    = params.mods;
smin    = params.smin;
smax    = params.smax;
gmin    = params.gmin;
gmax    = params.gmax;
N       = params.N;
walls   = params.walls;
dynamics= params.dynamics;
Nmods   = length(mods);


% ...
T = defaulter('sequencelength',T,varargin{:});
Ntraj = floor(Nexamples/T);
[Z,U,KFparams] = getCenterOutReacher(yrclass,Ntraj,T,dynamics,Ndims,mods,smin,smax);
keyboard
X = cat(2,Z,U);


%%%% CHECK/FIX ME
% store output and input for encoding in PPCs
%%% Y = state2stim(M,X,smin,smax,N,walls);
M = blkdiag(dynamics.C,dynamics.H);
Y = shortdata(Ntraj,4,(M*longdata(cat(2,Z,U))')');
Y = reshape(Y,[Ntraj,Ndims,Nmods,T]);
S(:,1:Ndims,strcmp(mods,'Hand-Position'),:) = Y(:,:,1,:);
S(:,1:Ndims,strcmp(mods,'Hand-Velocity'),:) = Y(:,:,2,:);
S(:,1:Ndims,strcmp(mods,'Efference-Copy'),:) = Y(:,:,3,:);
%%%%

S = state2stim(M,X,smin,smax,N,walls);

% output the stimuli, the state, locations of traj. restarts, and gains
S           = longdata(S);
Q.states    = longdata(X);
Q.restarts  = 1:T:(Ntraj*T);
Q.G = UniformNormalDiracSampler(zeros(size(gmin),yrclass),Inf,Ntraj*T,...
    gmin,gmax,0);
Q.KFparams = KFparams;



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function S = state2stim(M,X,smin,smax,N,walls)
% multiply by the emission matrix and then wrap

% params
[Ndims,Nmods] = size(smin);
[Ntraj,~,T] = size(X);
smin = shiftdim(smin,-1);
smax = shiftdim(smax,-1);

% turn into emissions
Sunwrapped = shortdata(Ntraj,4,reshape(longdata(X)*M',[Ntraj*T,Ndims,Nmods]));

% if requested, "wrap" the stimuli
if strcmp(walls,'wrapping')
    srange = N/(N-1)*(smax - smin);
    S = mod(Sunwrapped - smin,srange) + smin;
else
    S = Sunwrapped;
end

end
%-------------------------------------------------------------------------%