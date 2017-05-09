% Test if the rEFH uses the total spike counts to weight the "sensory"
% information.  
%
% Ideally, one would just run the normal test under differet, fixed gains.
% However, the total spike counts are too variable for any given gain; in
% particular, it is not uncommon to achieve higher total spike counts under
% the minimal gain than a non-trivial proportion of the total spike counts
% produced by the maximal gain.
%
% Therefore, ...

%-------------------------------------------------------------------------%
% Created: 01/05/14
%   -by JGM
%-------------------------------------------------------------------------%


% load
clear; clc;
load('dynamical\finalwts\wts1DrEFHManyXprmts.mat');
load('dynamical\finalwts\testdata1DrEFH.mat','Rtest','Xtest','Qtest');
LDSparamsTrue = getLDSparams(params.dynamics);

% decode inputs
[ShatTest,ttlSpksTest] = decodeDataPPC(Rtest,Xtest,Qtest,params);
InfoTest = GTPNposteriorInfo(ttlSpksTest,params);
pSENSORY = cumulantNeutralize(ShatTest,InfoTest,params);


% params
Ngains = 10;
Nruns = 12;
%%%%params.Ncases = 320; % there's too much variance in MSEs across tests
testgains = linspace(params.gmin,params.gmax,Ngains);

% malloc
statMat = zeros(Ngains,2,Nruns);

% loop
S = latents2stims(Xtest,Qtest.latent2stim,params.mods,params.Ndims);
NSind = strcmp(params.mods,params.NS);
for iGain = 1:Ngains
    for iRun = 1:Nruns
        
        params.universalGain = testgains(iGain);
        
        % filter with rEFH
        [~,~,Shat1,Info1] = testEFHPPC(Rtest,Xtest,Qtest,Allwts{1},params);
        pREFH = cumulantNeutralize(Shat1,Info1,params);
        err = pREFH.Xpct(:,:,NSind) - S(:,:,NSind);
        statMat(iGain,1,iRun) = err'*err/size(err,1);
        
        % Kalman filter with optimal parameters
        pOPT = KFposteriorization(pSENSORY,Q,LDSparamsTrue,'opt');
        err = pOPT.Xpct(:,:,NSind) - S(:,:,NSind);
        statMat(iGain,2,iRun) = err'*err/size(err,1);
    end
end

% store
xpctMat = mean(statMat,3);
stdMat = sqrt(var(statMat,[],3));

%%
%%% load results/sensitivityToInstantaneousReliability
f1 = figure; boxplot(squeeze(statMat(:,1,:))'); 
set(gca,'Xticklabel',sprintf('%1.2f\n ',testgains)); 
xlabel('gains'); ylabel('MSE (rad$^2$)');
matlab2tikzWrapper(['rEFHreliabilitySensitivity',date],f1);

f2 = figure; boxplot(squeeze(statMat(:,2,:))'); 
set(gca,'Xticklabel',sprintf('%1.2f\n ',testgains));
xlabel('gains'); ylabel('MSE (rad$^2$)');
matlab2tikzWrapper(['OPTreliabilitySensitivity',date],f2);













