function LDSparams = getLDSparams(dynamics)
% getLDSparams  Convert from params.dynamics to LDSparams
%
% USAGE:
%   LDSparams = getLDSparams(params.dynamics)
%
% This function does very little, mostly just copying over fields from one
% structure to another.  It also:
%   (1) puts the priors together into one field, and
%   (2) make sure a field muX exist (zeros if not a field in dynamics).


%-------------------------------------------------------------------------%
% Revised: 01/09/17
%   -reduced to just the case "true," as part of the Grand Revision
% Revised: 08/26/14
%   -changed case "true" to treat "efference copied" controls as states
% Revised: 01/10/14
%   -getXandY -> getLDSdata, and associated changes, esp. to accommodate
%   controlled systems
% Revised: 10/30/13
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%


%%% TO DO:
% (1) if there's a visible control....
% (2) GPUify based on params.machine??
% (3) get rid of muX0, muV0 in params.dynamics?
%%%%%%%%%%%%%%

% assemble prior
mu0     = [dynamics.muX0; dynamics.muV0];
Sigma0  = blkdiag(dynamics.SigmaX0,dynamics.SigmaV0);
if isfield(dynamics,'muU0')
    mu0     = [mu0; dynamics.muU0];
    Sigma0  = blkdiag(Sigma0, dynamics.SigmaU0);
end

% store
LDSparams.A = dynamics.A;
LDSparams.C = dynamics.C;
if isfield(dynamics,'muX')
    LDSparams.muX = dynamics.muX;
else
    LDSparams.muX = zeros(size(mu0));
end
LDSparams.SigmaX = dynamics.SigmaX;
LDSparams.mu0 = mu0;
LDSparams.Info0 = inv(Sigma0);
    
    
end