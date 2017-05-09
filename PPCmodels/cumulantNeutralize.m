function cumulants = cumulantNeutralize(xpctL,InfoL,params)
% cumulantNeutralize    Convert cumulants into "neutral space"
%
% USAGE:
%   posteriorCumulants = cumulantNeutralize(shatL,InfoL,params)
%
% For the underlying math, see:
%   [Makin, Fellows, Sabes; PLoS C.B., 2013], Supplementary Material
% The inputs 
%
%   xpctL (Nexamples x Ndims x Nmods) 
%
% and
%
%   InfoL  (Nexamples x Ndims x Ndims x Nmods) 
%
% are tensors containing the first two cumulants of a distribution.  But
% not all "mods" (third and fourth ways, resp.) are in the same space. This
% function converts all the cumulants to the "neutral space" given by
% params.NS, storing the results as fields in posteriorCumulants.
%
% This function is a countepart to getStimuliMultisensoryIntegration, and
% runs on basically the same datatypes.
%
% Note that for coordinate transformations (1Daddition, 2Daddition), this
% take in three modalities but outputs two, since two are combined.

%-------------------------------------------------------------------------%
% Revised: 01/04/17
%   -generalized to avoid case statements based on datatype; now parallels
%   getStimuliMultisensoryIntegration.m.
% Revised: 07/06/14
%   -added code for Gaussian prior
% Revised: 07/05/14
%   -added case '1Daddition'
% Created: 07/04/14
%   by JGM
%-------------------------------------------------------------------------%

% params
mods        = params.mods;
modNeutral  = params.NS;
if isfield(params,'roboparams') 
    roboparams = params.roboparams; 
else
    roboparams = [];
end

% a silly hack
div         = char(ones(1,length(mods))*'-');

% indices: "neutral," non-"neutral," and Gaze-Angle spaces
indNeutral  = strcmp(modNeutral,mods);
indGazeang  = strcmp('Gaze-Angle',mods);
indNonntrl  = ~indNeutral&~indGazeang;          % sum(.) >= 1
modNonntrl  = mods(indNonntrl);

% extract estimates (and info) from the respective modalities
shatNeutral = xpctL(:,:,indNeutral);
shatNonntrl = xpctL(:,:,indNonntrl);            % size(.,3) >= 1
shatGazeang = xpctL(:,:,indGazeang);
infoNonntrl = InfoL(:,:,:,indNonntrl);          % size(.,4) >= 1
infoGazeang = InfoL(:,:,:,indGazeang);
if ~any(indGazeang), shatGazeang = 0; end       % for '-Dinteg'

% store cumulants from the neutral modality
xpctNeutral(:,:,1)   = shatNeutral;
InfoNeutral(:,:,:,1) = InfoL(:,:,:,indNeutral);
srcs{1} = params.NS;

% store cumulants from the non-neutral modalities
for jMod = 1:sum(indNonntrl)
    
    % get the appropriate fxns for *this* non-neutral modality
    [Finv,Jfxn,InfoNeutralize] = setSensoryTransformations(modNeutral,...
        modNonntrl{jMod},roboparams);
    
    % store the cumulants
    xpctNeutral(:,:,1+jMod)   = Finv(shatNonntrl,shatGazeang);
    InfoNeutral(:,:,:,1+jMod) = InfoNeutralize(infoNonntrl(:,:,:,jMod),...
        infoGazeang,Jfxn(shatNeutral,shatGazeang));
    srcs{1+jMod} = [modNonntrl{jMod},div(indGazeang),mods{indGazeang}];
end
        

% add prior distribution, if there is one
if isfield(params,'p0')
    xpct0 = params.p0.mu;
    Info0 = inv(params.p0.cov);
    
    xpctNeutral(:,:,end+1) = xpct0(:,ones(size(xpctNeutral,1),1))';
    InfoNeutral(:,:,:,end+1) = permute(Info0(:,:,ones(size(xpctNeutral,1),1)),[3,1,2]);
    srcs{end+1} = 'prior';
end


% store
cumulants.Xpct = xpctNeutral;
cumulants.Info = InfoNeutral;
cumulants.srcs = srcs;


end


%%% You might consider allowing Gaze-Angle to be a "neutral space"--after
%%% all, that may be true in posterior parietal cortex!--but that would
%%% involve a bit more work.
% case 'Gaze-Angle'
% xpctNeutral(:,:,1) = ehat;
% xpctNeutral(:,:,2) = xhat - FK2link(thhat,roboparams,0);
%
% InfoNeutral(:,:,:,1) = eInfo;
% J = IKjacobianFast(xhat-ehat,roboparams);
% InfoNeutral(:,:,:,2) =...
% 1./(1./(J.*squeeze(thInfo).*J) + 1./squeeze(xInfo));

