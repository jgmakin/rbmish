function [LDSparamsEM, Bxz_EM, LDSparamsObs, XhatStatic, RsqStatic,...
    XhatDynamic, RsqDynamic] = latentStateLDSDecoder(getTrainingData,...
    Ytrain,Xtrain,Ytest,Xtest,Nmsperbin,NhiddenStates,LDSparamsEM,SStot)
% latentStateLDSDecoder
% USAGE: 
%   [LDSparamsEM, Bxz_EM, LDSparamsObs, XhatStatic, RsqStatic,...
%       XhatDynamic, RsqDynamic] = latentStateLDSDecoder(getTrainingData,...
%       Ytrain,Xtrain,Ytest,Xtest,Nmsperbin,NhiddenStates,LDSparamsEM,SStot);
%   

%-------------------------------------------------------------------------%
% Revised: 09/29/17
%   ......
% Cribbed: 09/29/17
%   -from filtersForNeuralData.m (JGM)
% Created: ~04/xx/16
%   by JGM
%-------------------------------------------------------------------------%

% learn a new system?
if isempty(LDSparamsEM)
    
    NepochsMax = 100;
%     %%%% Gharamani's code
%     tic
%     net = lds(gather(Ytrain),NhiddenStates,size(Ytrain,1),NepochsMax,0.0001);
%     toc
%     %%%%
    tic;
    LDSparamsEM = EM4LDS(NhiddenStates,NepochsMax,'Normal',...
        getTrainingData,...
        'verbosity', 1,...
        'diagonal covariances',[1,1,0],...
        'IC domain','all steps',...
        'parameter initialization','FactorAnalysis',...
        'convergence tolerance',0.0001);
    toc
end

% "observed variables"
LDSparamsEM.T   = size(Ytrain,1);
KFdstrbsTrainEM = KalmanFilter(LDSparamsEM,double(Ytrain'),'lightweight',1);
Vtrain          = [KFdstrbsTrainEM.XHATMU',ones(LDSparamsEM.T,1)];
clear KFdstrbsTrainEM Rtrain
LDSparamsEM.T   = size(Ytest,1);
KFdstrbsTestEM  = KalmanFilter(LDSparamsEM,double(Ytest'),'lightweight',1);
Vtest           = [KFdstrbsTestEM.XHATMU', ones(LDSparamsEM.T,1)];
clear KFdstrbsTestEM

% decode
[XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,Bxz_EM,LDSparamsObs] =...
    filterDecoder(Vtrain,Xtrain,Vtest,Xtest,1,SStot,Nmsperbin);

end
