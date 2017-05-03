function [thisDoB,stepsEMPAvg] = recalibrate(Di,IntegL0,wts,params)
% idea: The ostensible movement of the integrated estimated, from before to
% after recalibration, is actually the consequence of the model's visual
% space moving in the opposite direction.  So we can get an assay of the
% model's (changing) internal representation of visual space without any
% unimodal visual probe.

%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed the indexing of IntegL0 based on new indexing of shatL in
%   estStatsCorePP.m
% Revised: 09/27/12
%   -eliminated VWDR and WDR steps as output args, b/c they ought to be
%       calculated along their own paths rather than the empirical one;
%       this is now done in adaptByRule.m
% Revised: 09/19/12
%   -fixed some bugs
% Created: 09/18/12
%   by JGM
%   -this file is a heavily adapted version of recalibration2.m
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%
% NB: this function now expects Di,IntegL0 in longdata format and should be
% re-written accordingly (probably with no loops)
%%%%%%%%%%%%%%%%%%%


% init
Ndims = params.Ndims;
Nmods = length(params.mods);
[nCases,nDims,nBatches] = size(Di);


% malloc
stepsEMP    = zeros(Ndims,Nmods,nCases);
stepsEMPAvg = zeros(Ndims,Nmods,nBatches);

% outer loop (across "trials," i.e. training sessions)
tic
for iBatch = 1:nBatches
    
    thisDi = Di(:,:,iBatch);
    
    % integrate, retrain, then reintegrate
    thisDoA = updownDBN(thisDi,wts,params,'means','quiet'); % 'Nsamples','quiet');
    wts = train(thisDi,wts,params);
    thisDoB = updownDBN(thisDi,wts,params,'means','quiet'); % 'Nsamples','quiet');
    
    % get actual and predicted steps
    for iCase = 1:nCases
        
        % the data for this case (and this batch)
        di = thisDi(iCase,:);
        integL0 = IntegL0(iCase,:,:,iBatch);
        
        % *adjusted* unimodal ests, using prev wts
        shatNiA = getUnimodalEsts(di,thisDoA(iCase,:),integL0,params);
        
        % *adjusted* unimodal ests, using new wts
        shatNiB = getUnimodalEsts(di,thisDoB(iCase,:),integL0,params);
        
        % *actual* step taken!
        stepsEMP(:,:,iCase) = shatNiB - shatNiA;
        
    end
    
    % average the critical quantities so that you can plot them
    stepsEMPAvg(:,:,iBatch) = mean(stepsEMP,3);
    
end
toc



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function shatNi = getUnimodalEsts(di,do,integL0,params)

integL = decoder(do,params);                % nominal integrated estimate
cumAdpt = integL - integL0;                 % cumulative adaptation, local
shatLi = decoder(di,params) - cumAdpt;      % adjust by shifted rprstn
shatNi = estGather(shatLi,params);          % convert to neutral space

end
%-------------------------------------------------------------------------%