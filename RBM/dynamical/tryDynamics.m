%%% hard-coded for Ndims=1

clear; clc

USEEMMODELS = 0;


% set all params
params = setParams;
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;


% train (or just get)
LDSdata = getLDSdata(params);
LDSparamsObs = getLDSparams(params,'observed',LDSdata);
LDSparamsTrue = getLDSparams(params,'true');
if USEEMMODELS
    foo = params;
    %%% load('results/finalwts/wts1Dosc141009.mat')
    %%%
    load('C:\Users\makin\Desktop\caca\LDSparamsEM1DrEFHwithECRandInit2ndOrdXXX');
    LDSparamsEMBad = LDSparamsEM;
    load('C:\Users\makin\Desktop\caca\LDSparamsEM1DrEFHwithECRandInit3rdOrdXXX');
    LDSparamsEMGood = LDSparamsEM;
    %%%
    
    LDSparamsEMBad.T = params.dynamics.T;
    LDSparamsEMGood.T = params.dynamics.T;
    params = foo;
end

% test
LDSdataTest = getLDSdata(params);
[~,pSENSORY] = EFHfilter(LDSdataTest,[],params);
pKFtrue = KF4PPC(LDSdataTest,LDSparamsTrue,'opt');
pKFobs = KF4PPC(LDSdataTest,LDSparamsObs,'obs');
if USEEMMODELS
    pKFEMBad = KF4PPC(LDSdataTest,LDSparamsEMBad,'EM');
    pKFEMGood = KF4PPC(LDSdataTest,LDSparamsEMGood,'EM');
end


% train and test 
if isfield(params.dynamics,'F')
    LDSparamsObsNoCtrlDynamics = learnfullyobservedLDS(LDSdata,2,1,1,1);
    [pKFobsNoU, LDSparamsObsNoU] = obsKFwithNoInput(LDSdata,LDSdataTest,params);
    pKFobsNoCtrlDynamics = KF4PPC(LDSdataTest,LDSparamsObsNoCtrlDynamics,'obs');
end




% plot params
setColors;
dt = params.dynamics.dt;
T = params.dynamics.T;
tvec = 0:dt:(dt*T-dt);

% states
for iCase = 1:params.Ncases
%%%iCase = ceil(params.Ncases*rand);
    
figure(83); clf;        

% position
subplot(3,1,1); hold on;
title('Position')
plot(tvec,squeeze(LDSdataTest.Z(iCase,1,:))','b')
plot(tvec,squeeze(LDSdataTest.Y(iCase,1,:))','color',PROPcolor)
if isfield(params.dynamics,'F')
    plot(tvec,squeeze(pKFobsNoU.Xpct(iCase,1,:))','color',OBScolor) 
    plot(tvec,squeeze(pKFobsNoCtrlDynamics.Xpct(iCase,1,:))','color',EMcolor)
end
plot(tvec,squeeze(pKFtrue.Xpct(iCase,1,:))','color',OPTcolor)
if params.Ndims == 1
    plot(tvec,params.smin(1,strcmp(params.NS,params.mods))*ones(size(tvec)),'k');
    plot(tvec,params.smax(1,strcmp(params.NS,params.mods))*ones(size(tvec)),'k');
    plot(tvec,(params.smin(1,strcmp(params.NS,params.mods))+0.1)*ones(size(tvec)),'k:');
    plot(tvec,(params.smax(1,strcmp(params.NS,params.mods))-0.1)*ones(size(tvec)),'k:');
end
if USEEMMODELS
    plot(tvec,squeeze(pKFEMBad.Xpct(iCase,1,:))','color',VIScolor)
    plot(tvec,squeeze(pKFEMGood.Xpct(iCase,1,:))','color',EMcolor)
end
hold off;

pause(0.1);

end


% velocity
subplot(3,1,2); hold on;
title('Velocity')
plot(tvec,squeeze(LDSdataTest.Z(iCase,2,:))','b')

% control
if isfield(params.dynamics,'F')
    subplot(3,1,3); hold on;
    title('Control')
    plot(tvec,squeeze(LDSdataTest.Z(iCase,3,:))','b')
    plot(tvec,squeeze(LDSdataTest.Y(iCase,2,:))','color',PROPcolor)
    plot(tvec,squeeze(pKFobsNoU.Xpct(iCase,2,:))','color',OBScolor)
    plot(tvec,squeeze(pKFobsNoCtrlDynamics.Xpct(iCase,2,:))','color',OBScolor)
    plot(tvec,squeeze(pKFtrue.Xpct(iCase,2,:))','color',OPTcolor)
    hold off;
end

% the components that make up the velocity
figure(5782); clf; 
subplot(2,1,1); hold on;
title('components of position')
posVar = params.dynamics.SigmaX(1,1);
plot(tvec,squeeze(LDSdataTest.Z(iCase,1,:))','r');
plot(tvec,dt*squeeze(LDSdataTest.Z(iCase,2,:))','g');
plot(tvec,sqrt(posVar)*randn(size(tvec)),'m');
hold off;

% the components that make up the position
% subplot(2,1,2); hold on;
% title('components of velocity')
% a1 = LDSparamsTrue.A(2,2);
% a2 = LDSparamsTrue.A(3,3);
% m = 1/LDSparamsTrue.A(2,3);
% velVar = params.dynamics.SigmaX(2,2);
% plot(tvec,a1*squeeze(LDSdataTest.Z(iCase,2,:))','r');
% plot(tvec,squeeze(LDSdataTest.Z(iCase,3,:))'/m,'g');
% plot(tvec,sqrt(velVar)*randn(size(tvec)),'m');
% hold off;


% MSEs
if isfield(params.dynamics,'F')
   posteriors = {pSENSORY,pKFobsNoU,pKFobsNoCtrlDynamics,pKFobs,pKFtrue};
%     posteriors = {pSENSORY,pKFobsNoU,pKFEMBad,pKFEMGood,...
%         pKFobsNoCtrlDynamics,pKFobs,pKFtrue,};
else
    if USEEMMODELS
        posteriors = {pSENSORY,pKFobs,pKFtrue,pKFEMBad,pKFEMBad,pKFEMGood};
    else
        posteriors = {pSENSORY,pKFobs,pKFtrue};
    end
end
eStats = testDynamics(LDSdataTest,params,1,posteriors{:});
















