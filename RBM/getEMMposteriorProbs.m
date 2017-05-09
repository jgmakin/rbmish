function [Pizq,Thyq] = getEMMposteriorProbs(R,wts,params)
% Get posterior probs. and std. params under EFH Erlang Mixture Model
%
% USAGE:
%   [Pizq,Thyq] = getEMMposteriorProbs(R,wts,params)

%-------------------------------------------------------------------------%
% Revised: 12/07/16
%   -changed to take R in longdata, rather than shortdata, format
% Cribbed: 06/02/16
%   from evaluateEMMEFH.m
%   by JGM
%-------------------------------------------------------------------------%


% init
visDstrbs = params.typeUnits{1};
hidDstrbs = params.typeUnits{2};
visNums = params.numsUnits{1};
hidNums = params.numsUnits{2};
Nvis = visNums(end);

endinds = cumsum(hidNums);
startinds = [1, endinds(1:end-1)+1];

% get (conditional) std. params. for model hiddens
if strcmp(params.algorithm,'rEFH')
    testData.Y = 0;
    testData.SigmaYX = 0;
    %%%%%%% FIX ME!
    testData.restarts = cell(1,params.dynamics.T);
    testData.R = shortdata(size(R,1)/params.dynamics.T,3,R);
    %%%%%%%
    [~,Thzq] = EFHfilter(testData,wts,params);
    Thzq = longdata(Thzq);
else
    Thzq = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),...
        hidDstrbs,hidNums,params);
end

% get (conditional) probabilities for model hiddens
Pizq = [];
for iGrp = 1:length(hidDstrbs)
    thisThzq = Thzq(:,startinds(iGrp):endinds(iGrp));
    Pizq = cat(2,Pizq,getCatProbs(thisThzq,hidDstrbs{iGrp}));
end

% get (conditional) std. params for model visibles
Thhat = invParamMap(Thzq,wts{2}(1:end-1,:),wts{2}(end,:),...
    visDstrbs,visNums,params);
if strcmp(params.algorithm,'rEFH')
    Nrcrnt = params.numsUnits{1}(1);
    if strcmp(params.typeUnits{1}(end),'GammaFixedScale')
        Thyq(:,1) = mean(Thhat(:,Nrcrnt+(1:Nvis)),2);
        Thyq(:,2) = params.scaleparams;
    else 
        Thyq(:,1) = mean(Thhat(:,Nrcrnt+(1:Nvis/2)),2);
        Thyq(:,2) = mean(Thhat(:,Nrcrnt+((Nvis/2+1):Nvis)),2);
    end
else
    if strcmp(params.typeUnits{1}(end),'GammaFixedScale')
        Thyq(:,1) = mean(Thhat(:,1:Nvis),2);
        Thyq(:,2) = params.scaleparams;
    else
        Thyq(:,1) = mean(Thhat(:,1:Nvis/2),2);
        Thyq(:,2) = mean(Thhat(:,(Nvis/2+1):Nvis),2);
    end
end


end






