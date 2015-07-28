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
load('results/finalwts/wts1DrEFH141029.mat')
%%% comment this out in LDStest.m !!! as well as clear; clc;

% params
Ngains = 10;
Nruns = 12;
params.Ncases = 320; % there's too much variance in MSEs across tests
gmin = (1 - params.swing)*params.g;
gmax = (1 + params.swing)*params.g;
testgains = linspace(gmin,gmax,Ngains);

% malloc
statMat = zeros(Ngains,2,Nruns);

% loop
for iGain = 1:Ngains
    for iRun = 1:Nruns
        
        params.universalGain = testgains(iGain);
        
        LDStest
        posteriorNames = {eStats.tags(:).name};
        inds = logical(strcmp(posteriorNames,'rEFH') +...
            strcmp(posteriorNames,'opt'));
        statMat(iGain,:,iRun) =...
            eStats.Xpct(inds).^2 + squeeze(eStats.Cvrn(inds))';
        
    
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













