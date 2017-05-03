clear; clc; close all;
load adpt0315avg

nBatches = 120;
stddevs = 5;
[Di,Si,shft] = generatebiaseddata(nBatches,stddevs,params);

[stVWDR,stWDR] = adaptByRule(Di,squeeze(Si(1,:,:)),50,0.01,1000,params);
plotAdaptations(Si,IntegL0(:,:,end),DoBF,stVWDR,stWDR,stEMP,params);

figure; hold on
plot(squeeze(stVWDR(2,2,:)))
plot(squeeze(stWDR(2,2,:)),'r')
hold off

% save adpt0315avg IntegL0 DoBF stVWDR stWDR stEMP params