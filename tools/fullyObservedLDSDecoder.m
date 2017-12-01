function [LDSparamsOBS, Xhat_kf, RsqKF, Xhat_rtss, RsqRTSS] =...
    fullyObservedLDSDecoder(Ytrain,Xtrain,T,Ytest,Xtest,SStot)
% fullyObservedLDSDecoder
%
% USAGE:
% [LDSparamsObs, KFdstrbsObs, RTSSdstrbsObs,...
%   Xhat_kf, RsqKF_Obs, Xhat_rtss, RsqRTSS_Obs] = ...
% 	fullyObservedLDSDecoder(Rtrain,Xtrain,Qtrain.T,Rtest,Xtest,SStot);

%-------------------------------------------------------------------------%
% Cribbed: 09/29/17
%   -from filtersForNeuralData.m (JGM)
% Created: ~04/xx/16
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nexamples = size(Ytrain,1);
Ntraj = floor(Nexamples/T);

% train
LDSparamsOBS = learnfullyobservedLDS(shortdata(Ntraj,3,Ytrain),...
    shortdata(Ntraj,3,Xtrain));

% but use *all* data to get initial state (overwrite)
LDSparamsOBS.mu0 = mean(Xtrain)';
LDSparamsOBS.Info0 = inv(cov(Xtrain));

% test: filter, smoother
tic;
KFdstrbs    = KalmanFilter(LDSparamsOBS,gather(Ytest'));
Xhat_kf     = KFdstrbs.XHATMU';
SSerr       = sum((Xhat_kf - Xtest).^2);
RsqKF       = 1 - SSerr./SStot;
RTSSdstrbs  = RTSsmoother(LDSparamsOBS,KFdstrbs);
Xhat_rtss   = RTSSdstrbs.XHAT';
SSerr       = sum((Xhat_rtss - Xtest).^2);
RsqRTSS     = 1 - SSerr./SStot;
toc;


end
