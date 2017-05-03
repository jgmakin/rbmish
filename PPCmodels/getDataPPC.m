function [R,Q] = getDataPPC(X,Q,params,varargin)
% getDataPPC    Linear probabilistic population code response to stimuli
% 
% USAGE:
%   R = getDataPPC(X,Q,params)
%
% Given latent variables X and the data structure Q, getDataPPC generates 
% training/testing data R for EFHs.  The data are "neural responses" to the
% "stimuli" (possibly transformed versions of X) for populations of
% Gaussian-tuned neurons that smoothly tile their respective spaces.  The 
% popluation gains should be a field G in Q.  Other parameters must be in 
% the params structure.  For inputs of size:
%
%       S: (Nexamples x Ndims x Nmods)
%       G: (Nexamples x Nmods)
%
% the function returns:
%
%       R: (Nexamples x Nmods*N^Ndims)
 

%-------------------------------------------------------------------------%
% Cribbed: 12/31/16
%   by JGM
%   from generateData.m
%-------------------------------------------------------------------------%

% params
visDstrbs   = params.typeUnits{1};
visNums     = params.numsUnits{1};
mods        = params.mods;
walls       = params.walls;
Ndims       = params.Ndims;


% set the biases (default to none)
for iMod = 1:length(mods)
    jMod = strcmp(mods,mods(iMod));     % *could* have >1 nonzero entry
    biases(:,jMod) = defaulter(['bias: ',mods{iMod}],...
        zeros(size(params.smin,sum(jMod)),1,'like',X),varargin{:});
end

% convert latent state to "stimuli" to be encoded
S = latents2stims(X,Q.latent2stim,mods,Ndims);

if strcmp(walls,'wrapping')
    R = encodeToroidalStimuli(S,Q.G,biases,visDstrbs,visNums,params);
else
    R = stdStimulusEncoding(S,Q.G,biases,visDstrbs,visNums,params);
end

% mark the decoupled vectors
%%%R(:,params.N) = Q.DECOUPLED*params.gmax(1) + ~Q.DECOUPLED.*R(:,params.N);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function R = encodeToroidalStimuli(S,G,biases,dstrbs,nums,params)

% wrap the stimuli into encoding space
[S,srange] = wrapStimuli(S,params.smin,params.smax,params.N);

% then displace them by length of range(s) for the recursive thing
S = S - srange;

% since all the gains have been sampled, don't resample
[R,~] = encodeRecursively(S,G,biases,dstrbs,nums,1,params);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [R,S] = encodeRecursively(S,G,biases,dstrbs,nums,dim,params)

% init
R = 0;
[~,Ndims,Nmods] = size(S);
N = params.N;

% (re-)calculate srangemat
srange = reshape(N/(N-1)*(params.smax - params.smin),[1,Ndims,Nmods]);

if dim <= Ndims
    thisS = S;
    for iWrap = 1:3
        [thisR,thisS] = encodeRecursively(thisS,G,biases,dstrbs,nums,dim+1,params);
        R = R + thisR;
        
        thisS(:,dim,:) = thisS(:,dim,:) + srange(1,dim,:);
    end
else
    R = stdStimulusEncoding(S,G,biases,dstrbs,nums,params);
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function R = stdStimulusEncoding(S,G,biases,dstrbs,nums,params)

% Ns
[Nexamples,Ndims,Nmods] = size(S);
N = params.N;

% draw population responses
R = zeros(Nexamples,Nmods*N^Ndims,'like',S);
for iMod = 1:Nmods
    
    % prepare
    b = biases(:,iMod);
    inds = ((iMod-1)*N^Ndims + 1):(iMod*N^Ndims);
   
    % encode the stimuli in GTPNs
    R(:,inds) = PPCencode(S(:,:,iMod) + b',G(:,iMod),...
        params.smin(:,iMod),params.smax(:,iMod),dstrbs,nums,...
        params);
    
end


end
%-------------------------------------------------------------------------%