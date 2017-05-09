function F = tiledTuning(S,g,smin,smax,dstrbs,nums,params)
% tiledTuning   A population of smoothly tiled tuning curves
%
% USAGE:
%   F = tiledTuning(X,g,smin,smax,inputDstrb,params);
%
% Produces a matrix F of size (Nexamples x Nneurons) when given a matrix of
% inputs S (Nexamples x Ndims), a vector of gains g (Nexamples x 1), the
% range of the stimuli, smin and smax, and the params structure.
%
% NB: main use is in PPCencode, but it can be used on its own for a
% "noiseless responses."

%-------------------------------------------------------------------------%
% Revised: 01/02/17
%   -added nums argument and made dstrbs a cell array of strings rather
%   than just a string.  Now this function has to loop through those
%   entries, which is "more general."
% Revised: 08/23/16
%   -moved rescaling from generateData into this function
% Revised: 08/22/16
%   -added linear tuning (no preferred directions)
% Revised: 02/23/16
%   -changed to accommodate *array* of inputUnitTypes.
%   -added numsUnits as an argument
% Revised: 08/07/14
%   -replaced machine ('domestica') check with data class check
% Revised: 07/29/14
%   -added gpu stuff
% Revised: 07/24/14
%   -rewrote from scratch to use binary singleton expansion, tensorOp.m
%   -incorporated Ndims = 1, 2 to one case
%   -generalized to any Ndims
% Adapted: 06/21/14
%   -from GTrespfxn (version history below)
% Revised: 05/06/14
%   -vectorized inputs and outputs! (avoids parfor loop)
%   -renamed from respfxn to GTrespfxn
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) make the tuning curves a seperate parameter??
% (2) This technically accepts cell arrays dstrbs that contain more than
% one kind of distribution--but then it just uses the same gains g and the
% same smin and smax for all of these.  That doesn't seem useful....


if onlyOneUniqueStr(dstrbs)
    F = stimuli2tuningCurves(S,g,smin,smax,dstrbs{1},params);
else
    endinds = cumsum(nums);
    startinds = [1, endinds(1:end-1)+1];
    for iGrp = 1:length(dstrbs)
        F(:,startinds(iGrp):endinds(iGrp)) = stimuli2tuningCurves(...
            S(:,startinds(iGrp):endinds(iGrp)),g,smax,smin,dstrbs{iGrp},...
            params);
    end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function F = stimuli2tuningCurves(S,g,smin,smax,dstrb,params) 

switch dstrb
    case 'StandardNormal'   % linear tuning
        F = g.*S;
        %%% NB that this ignores all tuning structure in params
        
    case {'Bernoulli','Binomial','BernoulliDropout'}
        
        % rescale so that stimuli span [0,1] x [0,1] x ...
        patchmin = params.margin*ones(1,params.Ndims);
        patchmax = patchmin + params.respLength;
        X = scalefxn(S,smin,smax,patchmin,patchmax);
        
        % subtract from logit(g), send through the logistic function
        logit = @(z)(log(z./(1-z)));
        logistic = @(z)(1./(1 + exp(-z)));
        quadraticForm = quadraticTune(X,params.C,params.N,params.gridsize);
        F = logistic(logit(g) - quadraticForm/2);
        
    otherwise % Gaussian tuning
        
        % rescale so that stimuli span [0,1] x [0,1] x ...
        patchmin = params.margin*ones(1,params.Ndims);
        patchmax = patchmin + params.respLength;
        X = scalefxn(S,smin,smax,patchmin,patchmax);
        %%% what's with the transpose??
        
        % exponentiate quadratic form to get a Gaussian, scale by gains
        quadraticForm = quadraticTune(X,params.C,params.N,params.gridsize);
        F = g.*exp(-quadraticForm/2);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function quadraticForm = quadraticTune(X,C,N,gridsize)
% compute (x - pd)'*Info*(x - pd)

% params
Info = inv(C);
[Nexamples, Ndims] = size(X);

% build a lattice across Ndims with points at the PDs
PDs = linspace(0,gridsize,N)*ones(1,1,'like',X);
latticePDs = ndlattice(Ndims,PDs);

Infotensor = repmat(cast(Info,'like',X),[1,1,Nexamples]);
Y = shiftdim(X',-1) - latticePDs;           % Nunits x Ndims x Nexamples
YtrS = tensorOp(Y,Infotensor);              % Nunits x Ndims x Nexamples
quadraticForm = permute(sum(YtrS.*Y,2),[3,1,2]);    % Nexamples x Nunits

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function FLAG = onlyOneUniqueStr(strcell)

d = strcell{1};
FLAG = true;
i=1;
while FLAG&&(i<length(strcell))
    FLAG = strcmp(strcell{i+1},d);
    i=i+1;
end

end
%-------------------------------------------------------------------------%
