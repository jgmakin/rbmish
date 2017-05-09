function [X,Q] = getLatentsMultisensoryIntegration(Nexamples,yrclass,...
    source,params,varargin)
% getLatentsMultisensoryIntegration
%
% USAGE:
%   [X,Q] = getLatentsMultisensoryIntegration(Nexamples,yrclass,...
%       params.NS,params);
%
% This function expects two or three entries in mods.  The tensor X has
% size (Nexamples x Nstates x Nltnts), where Nstates = Ndims for each of
% the Nltnts=2, the neutral space and the gaze angle (NB that we assume 
% 'Gaze-Angle' will never itself be the neutral space).  (A distinction is
% made between Nstates and Ndims for parallelism with getLatentsLTI.m).  To
% convert from X to *stimuli*--e.g., for encoding in PPCs, or for computing
% decoding errors--use the conversion functions stored in the cell array
% Q.latent2stim.

%-------------------------------------------------------------------------%
% Revised: 01/05/17
%   -Grand Revision still in progress
%   -changed outputs from stimuli to latent variables, to parallel 
% Revised: 12/30/16
%   -rewrote to run on inputs rather than fields in params and Q.
% Cribbed: 12/30/16
%   from generateData.m
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) Restore decoupling, which you removed.  The idea would be to make a
% third latent variable, after sourcevars and gaze angles, just in case 
% any(DECOUPLED).  Then the latent2[] fxns would have to change, too.


% params
Ndims   = params.Ndims;
mods    = params.mods;
smin    = params.smin;
smax    = params.smax;
gmin    = params.gmin;
gmax    = params.gmax;
if isfield(params,'roboparams')
    roboparams = params.roboparams;
else
    roboparams = [];
end

% indices
srcInd  = strcmp(source,mods);
eyeInd  = strcmp('Gaze-Angle',mods);

% set defaults
pX0.mu  = NaN(Ndims,1,yrclass);
pX0.cov = Inf; % => uniform distribution
pX0     = defaulter('latentvarprior',pX0,varargin{:});
crln    = defaulter('correlation',1,varargin{:});
pE0.mu  = zeros(Ndims,1,yrclass);
if any(eyeInd), pE0.cov = Inf; else pE0.cov = 0; end
pE0     = defaulter('gazeangleprior',pE0,varargin{:});
pG0.mu  = NaN(length(gmin),1,yrclass);
pG0.cov = Inf; % => uniform distribution

% speak!
describeData(Nexamples,Ndims,mods,crln,gmin,gmax,params.datatype);

% always draw samples of the source modality *and* gaze angles
sourcevars = defaulter('stims',UniformNormalDiracSampler(pX0.mu,pX0.cov,...
    Nexamples,smin(:,srcInd),smax(:,srcInd),0),varargin{:});
gazeangles = defaulter('gazes',UniformNormalDiracSampler(pE0.mu,pE0.cov,...
    Nexamples,smin(:,eyeInd),smax(:,eyeInd),0),varargin{:});
X = cat(3,sourcevars,gazeangles);

% build the functions that convert from latent vars to "stimuli"
Q.latent2stim = setConversionFxns(mods,source,roboparams);

% draw the gains
Q.G = defaulter('gains',UniformNormalDiracSampler(pG0.mu,pG0.cov,...
    Nexamples,gmin,gmax,0),varargin{:});


% 
%%% Q.DECOUPLED = rand(Nexamples,1) >= crln;

% check for out-of-range values of the stimulus
%%%%% put this back in!  It's just that, on the torus, you think you
%%%%% can get values slightly over the max......
%     for iMod = 1:size(S,3)
%         for iExample = 1:size(S,1)
%             s = squeeze(S(iExample,:,iMod));
%             outofrange(s,params.smin(:,iMod),params.smax(:,iMod),...
%                 params.mods{iMod},1,params);
%         end
%     end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function latent2stim = setConversionFxns(mods,source,roboparams)

Nmods = length(mods);
latent2stim = cell(1,Nmods);
for iMod = 1:Nmods
    f = setSensoryTransformations(mods{iMod},source,roboparams);
    latent2stim{iMod} = @(X)(f(X(:,:,1),X(:,:,2)));
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function describeData(Nexamples,Ndims,mods,crln,gmin,gmax,datatype)

fprintf('\n\n');
fprintf('Generating %i stimuli of type %s\n',Nexamples,datatype);
fprintf('\t%iD data in %i modalities\n',Ndims,length(mods));
fprintf('\tmodality correlation: %1.2f\n',crln);
for iMod = 1:length(mods)
    fprintf('\ttuning curve max between: [%1.2f,%1.2f]\n',gmin(iMod),gmax(iMod));
end
fprintf('\n\n');

end
%-------------------------------------------------------------------------%