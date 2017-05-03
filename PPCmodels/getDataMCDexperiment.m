function [R,Q] = getDataMCDexperiment(S,Q,params)
% getDataMCDexperiment    Linear probabilistic population code for MCD exp.
% 
% USAGE:
%   R = getDataMCDexperiment(S,Q,params)
%

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from generateData
%   by JGM
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




% pull out the stimuli
Sp = S(:,:,strcmp(params.mods,'ICMSpolar'));
Sc = S(:,:,strcmp(params.mods,'Motion-Dots'));

% encode the ICMS via electrodeFunc
F = params.experiment.electrodeFunc(Sp(:,1),Sp(:,2));
Nns = params.experiment.neuronsPerElectrode;
Ricms = sampleT(repmat(F,[1,1,Nns]),visDstrbs,visNums,params);
Ricms = reshape(Ricms,[size(F,1),size(F,2)*Nns]);
%%%%
% The icms gain is never used---that looks like a bug
%%%%


% encode the visual stimuli in the standard way
visInd = strcmp(params.mods,'Motion-Dots');
Rvis = PPCencode(Sc + Q.biases(:,visInd)',Q.G(:,visInd),...
    params.smin(:,visInd),params.smax(:,visInd),visDstrbs,visNums,params);

% store
R(:,:,strcmp(params.mods,'Motion-Dots'))    = Rvis;
R(:,:,strcmp(params.mods,'ICMSpolar'))      = Ricms;


end