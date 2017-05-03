function [Shat,ttlSpks,err] = decodeDataPPC(R,X,Q,params)
% GTPNdecode
%
% USAGE:
%   [Shat,ttlSpks,err] = decodeDataPPC(R,X,Q,params)

%-------------------------------------------------------------------------%
% Created: 01/07/17
%   by JGM
%-------------------------------------------------------------------------%

% only decode the Poisson parts of the visible layer; this only has an
% effect on certain models (e.g., 'HierL2').
visDstrbs   = params.typeUnits{1};
visNums     = params.numsUnits{1};
endinds     = cumsum(visNums);
startinds   = [1, endinds(1:end-1)+1];
iPPC        = strcmp(visDstrbs,'Poisson');
PPCinds     = startinds(:,iPPC):endinds(:,iPPC);
PPCs        = R(:,PPCinds);

% minimize sufficient statistics, compute errors
[Shat, ttlSpks] = GTPNsuffstats(PPCs,params);
Snotwrapped = latents2stims(X,Q.latent2stim,params.mods,params.Ndims);
if strcmp(params.walls,'wrapping')
    [Shat, err] = getBestTorusEsts(Shat,Snotwrapped,params);
else
    err = Shat - Snotwrapped;
end

end