function [stepVWDR,stepWDR,inputDscrpN] = getPredictedStep(di,tuningCov,...
    bias,shatLiA,shatNiA,params)

% input discrepancy: (v - p), *in prop space*, prior to adpt
inputDscrpN = shatNiA(:,1) - shatNiA(:,2);

% local covariances in prop space
SILSCerr = PPCinputStats(di,tuningCov,bias);
[~,SINSCondCov] = SICE(shatLiA(:)',params,SILSCerr);
covVN = squeeze(SINSCondCov(:,:,:,1));       %%% hard-coded "V"
covPN = squeeze(SINSCondCov(:,:,:,2));       %%% hard-coded "P"

% predictions according to VWDR and WDR
stepVWDR = [-covVN*inputDscrpN, covPN*inputDscrpN];
stepWDR  = [-covVN/(covVN + covPN)*inputDscrpN,...
    covPN/(covVN + covPN)*inputDscrpN];

end