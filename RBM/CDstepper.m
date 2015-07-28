function [visstates, hidstates] = CDstepper(hidstates,vishid,visbiases,...
    hidbiases,HIDFXN,VISFXN,params)
% contrastive divergence Gibbs sampler for RBMs
%-------------------------------------------------------------------------%
% Revised: 10/25/10
%   -changed to use all samples (no more "means")
% Adapted: 08/03/10 (happy b'day, TAM)
%   by JGM
%-------------------------------------------------------------------------%

% Gibbs sample
for step = 1:params.numCDsteps
    
    % feed down and back up
    vismeans = feedforward(hidstates,vishid',visbiases',VISFXN,params);
    visstates = sampler(vismeans,VISFXN,params);
    hidmeans = feedforward(visstates,vishid,hidbiases,HIDFXN,params);
    hidstates = sampler(hidmeans,HIDFXN,params);
    
    if sum(sum(isnan(vismeans)))
        fprintf('warning: confabulations contain NaNs\n');
        pause()
    end
   
end

end
