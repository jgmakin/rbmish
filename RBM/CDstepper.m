function [visstates, hidstates] = CDstepper(hidstates,vishid,visbiases,...
    hidbiases,hidDstrbs,visDstrbs,hidNums,visNums,Ncdsteps,params)
% contrastive divergence Gibbs sampler for RBMs
%-------------------------------------------------------------------------%
% Revised: 10/25/10
%   -changed to use all samples (no more "means")
% Adapted: 08/03/10 (happy b'day, TAM)
%   by JGM
%-------------------------------------------------------------------------%

% Gibbs sample
for step = 1:Ncdsteps
    
    % feed down and back up
    vismeans = invParamMap(hidstates,vishid',visbiases',visDstrbs,visNums,params);
    visstates = sampleT(vismeans,visDstrbs,visNums,params);
    hidmeans = invParamMap(visstates,vishid,hidbiases,hidDstrbs,hidNums,params);
    hidstates = sampleT(hidmeans,hidDstrbs,hidNums,params);
   
    if any(isnan(vismeans(:)))
        fprintf('warning: confabulations contain NaNs\n');
        pause()
    end
    
%     topo = displayshape(vismeans(2,:),params);
%     imagesc(topo{1})
%     pause;
   
end
%%%%%%
hidstates = hidmeans;
%%%%%%

end
