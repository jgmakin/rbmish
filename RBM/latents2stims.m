function S = latents2stims(X,latent2stimFxns,mods,Ndims)
% latents2stims     Convert latent variables into stimuli (to be encoded)
%
% USAGE:
%   S = latents2stims(X,Q.latent2stim,params.mods,params.Ndims)
%
% This functions is kind of like the "body": all the transduction of the
% world's (latent) variables prior to the neural responses

%-------------------------------------------------------------------------%
% Created: 01/05/17
%   -part of the Grand Revision 
%   -by JGM
%-------------------------------------------------------------------------%

% malloc
Nmods = length(mods);
S = zeros([size(X,1),Ndims,Nmods],'like',X);

% "stimuli"
for iMod = 1:Nmods
    S(:,:,iMod) = latent2stimFxns{iMod}(X);
end

end
