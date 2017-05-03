clear; clc; close all;

%% some spikes
setColors;
Nsamples = 20;
lambdaMax = 15;
mu = -2:3;

lambda = lambdaMax*exp(-mu.^2/2);
r = sampleT(repmat(lambda,Nsamples,1),'Poisson',length(lambda),[]);

hh = figure();
hold on
for iSample = 1:Nsamples
    scatter(1:6,r(iSample,:),[],repmat(OBScolor,6,1),'filled')
end
hold off;

ylabel('spikes')
set(gca,'XTick',[1:12]);

matlab2tikzWrapper(['someSpikes',date],hh);

%% multisensory integration schematic
close all;

x = linspace(-10,10,1000);
s1 = 1; mu1 = 5; Z1 = sqrt(2*pi*s1);
cue1 = exp(-(mu1 - x).^2/s1/2)/Z1;

s2 = 3; mu2 = -3; Z2 = sqrt(2*pi*s2);
cue2 = exp(-(mu2 - x).^2/s2/2)/Z2;

s3 = inv(inv(s1) + inv(s2));
mu3 = s3*inv(s1)*mu1 + s3*inv(s2)*mu2;
Z3 = sqrt(2*pi*s3);
integ = exp(-(mu3 - x).^2/s3/2)/Z3;

hh = figure; hold on;
plot(x,cue1,'color',VIScolor);
plot(x,cue2,'color',PROPcolor);
plot(x,integ,'color',OPTcolor);
hold off;

set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
ylabel('likelihood');
xlabel('$\Stim$','Interpreter','none');

matlab2tikzWrapper(['MIconceptual',date],hh);



%% hierarchical integration schematic
close all;

x = linspace(-10,10,1000);
s1 = 1; mu1 = 5; Z1 = sqrt(2*pi*s1);
cue1 = exp(-(mu1 - x).^2/s1/2)/Z1;

s2 = 3; mu2 = -3; Z2 = sqrt(2*pi*s2);
cue2 = exp(-(mu2 - x).^2/s2/2)/Z2;

s3 = 3; mu3 = 1; Z3 = sqrt(2*pi*s3);
cue3 = exp(-(mu3 - x).^2/s3/2)/Z3;


s4 = inv(inv(s1) + inv(s2) + inv(s3));
mu4 = s4*inv(s1)*mu1 + s4*inv(s2)*mu2 + s4*inv(s3)*mu3;
Z4 = sqrt(2*pi*s4);
integ = exp(-(mu4 - x).^2/s4/2)/Z4;

hh = figure; hold on;
plot(x,cue1,'color',VIScolor);
plot(x,cue2,'color',PROPcolor);
plot(x,cue3,'c');
plot(x,integ,'color',OPTcolor);
hold off;

set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
ylabel('likelihood');
xlabel('$\Stim$','Interpreter','none');

matlab2tikzWrapper(['conceptualHierL2',date],hh);


%% coordinate transformation schematic
close all;

x = linspace(-10,10,1000);
s1 = 1; mu1 = 5; Z1 = sqrt(2*pi*s1);
cue1 = exp(-(mu1 - x).^2/s1/2)/Z1;

s2 = 2; mu2 = -4.5; Z2 = sqrt(2*pi*s2);
cue2 = exp(-(mu2 - x).^2/s2/2)/Z2;

s3 = 3; mu3 = 1.5; Z3 = sqrt(2*pi*s3);
cue3 = exp(-(mu3 - x).^2/s3/2)/Z3;


s4 = inv(inv(s1+s2) + inv(s3));
mu4 = s4*inv(s1+s2)*(mu1+mu2) + s4*inv(s3)*mu3;
Z4 = sqrt(2*pi*s4);
integ = exp(-(mu4 - x).^2/s4/2)/Z4;

hh = figure; hold on;
plot(x,cue1,'color',VIScolor);
plot(x,cue2,'color',EYEcolor);
plot(x,cue3,'color',PROPcolor);
plot(x,integ,'color',OPTcolor);
hold off;

set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
ylabel('likelihood');
xlabel('$\Stim$','Interpreter','none');

matlab2tikzWrapper(['conceptual1Daddition',date],hh);


%% example trajectories
clear; clc; 
%%%load([getdir('data'),'RBMish/EFHs/wts_1D-LTI-PPC_rEFH_ManyXprmts'])
load([getdir('data'),'RBMish/EFHs/new/wts_rEFH_1D_LTI-PPC_240.mat'])
setColors;

if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
Ntraj = 1;
T = 1000;


% get data
[X,Q] = params.getLatents(Ntraj*T,dataclass);
R0 = params.getData(X,Q);

% decode inputs
[Shat0,ttlSpks0] = decodeDataPPC(R0,X,Q,params);
Info0 = GTPNposteriorInfo(ttlSpks0,params);
pSENSORY = cumulantNeutralize(Shat0,Info0,params);

% filter with rEFH
R1 = updownRDBN(R0,wts,params,T);
[Shat1,ttlSpks1] = decodeDataPPC(R1,X,Q,params);
Info1 = GTPNposteriorInfo(ttlSpks1,params); 
pREFH = cumulantNeutralize(Shat1,Info1,params);

% Kalman filter with true LTI sys
LDSparamsTrue = getLDSparams(params.dynamics);
pOPT = KFposteriorization(pSENSORY,Q,LDSparamsTrue,params); 


sampleRange = 551:720;
timeRange = sampleRange*params.dynamics.dt;

hh = figure(222); clf; hold on;
set(hh, 'Position',[100 493 500 300])
plot(timeRange,squeeze(X(sampleRange,1)),'color',[0.5,0.5,0.5],'linewidth',1.5);
plot(timeRange,squeeze(pSENSORY.Xpct(sampleRange,1)),'color',PROPcolor,'linewidth',1.5);
plot(timeRange,squeeze(pOPT.Xpct(sampleRange,1)),'color',OPTcolor,'linewidth',1.5);
plot(timeRange,squeeze(pREFH.Xpct(sampleRange,1)),'color',EFHcolor,'linewidth',1.5);
xlabel('time (s)');
ylabel('position (rad)');
hold off


matlab2tikzWrapper(['exampleTrajectories',date],hh);

%%
% clear; clc; close all
% load('results\finalwts\1DwrappingHVNdamped.mat')

setColors;

iTraj= 2;

hh = figure(101); clf;
time = linspace(0,5,1000);
plot(time,squeeze(LDSdataTest.Z(iTraj,1,:)),'color',[0.5 0.5 0.5])
xlabel('time (s)');
ylabel('position (rad)');
% matlab2tikzWrapper(['sample1dTraj',date],hh);


hold on;
plot(time,squeeze(pEFH.Shat(iTraj,1,1,:)),'color',EFHcolor);
hold off;


hold on;
plot(time,squeeze(pSENSORY.Shat(iTraj,1,1,:)),'color',PROPcolor);
hold off;


hold on;
plot(time,squeeze(pKFtrue.Shat(iTraj,1,1,:)),'color',OPTcolor);
hold off;

pause();

% end


%% reproduce MCD plots with one of your models
clear; clc
load results/finalwts/wts2Dallgains140725.mat

% params
testgains = [0.5,3,6.8,15,19];
%%% correspond to (scaled-versions of) the error variance in MCD plot

% test at different gains
errVrncEFH = zeros(length(testgains),4);
for iGain = 1:length(testgains)
    
    params.gmin = [testgains(iGain),20];
    params.gmax = [testgains(iGain),20];
    eStats = mastertest(wts,params);
    errVrncEFH(iGain,:) = arrayfun(@(i)(sqrt(det(eStats.Cvrn(:,:,i)))),...
        1:size(eStats.Cvrn,3));
end

% plot the standard deviations, a la MCD
fighandle = figure(4313); clf; hold on;
for iStat = 1:length(eStats.tags)
    shadedErrorBar(testgains,sqrt(errVrncEFH(:,iStat)),...
         0.71*sqrt(errVrncEFH(:,iStat))/sqrt(40000),...
         {'Color',getColor(eStats.tags(iStat).name)},2);
    %%% standard error of the standard deviation
    % plothandle = plot(testgains,sqrt(errVrncEFH(:,iStat)));
    % set(plothandle,'color',getColor(eStats.tags(iStat).name),'LineWidth',2);    
end
hold off;
xlabel('visual gain');  ylabel('std (cm)');
xticklabels = testgains;
xticks = testgains;
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
matlab2tikzWrapper(['2DallgainsMCDplot',date],fighandle);



% plot the standard deviations for MCD's data (Monkey D)
load EFH4MCD\iavalDmitri.mat
errStdD = iaval.unbiasederrstd.mean;
errStdStdD = iaval.unbiasederrstd.std;
errVrncEXP(:,1) = errStdD(:,1).^2;  errStdStdEXP(:,1) = errStdStdD(:,1);
errVrncEXP(:,2) = errStdD(:,3).^2;  errStdStdEXP(:,2) = errStdStdD(:,3);
errVrncEXP(:,4) = errStdD(:,2).^2;  errStdStdEXP(:,4) = errStdStdD(:,2);
errVrncEXP(:,3) = 1./(1./errVrncEXP(:,1) + 1./errVrncEXP(:,2));
errStdStdEXP(:,3) = 0.0001*ones(size(errStdStdD(:,3)));
coherences = [0,15,25,50,100];

fighandle = figure(4314); clf; hold on;
for iStat = 1:length(eStats.tags)
    shadedErrorBar(coherences,sqrt(errVrncEXP(:,iStat)),...
         errStdStdEXP(:,iStat),{'Color',getColor(eStats.tags(iStat).name)},2);
end
hold off;
xlabel('visual coherence');  ylabel('std (rad)');
xticklabels = coherences;
xticks = coherences;
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
matlab2tikzWrapper(['DmitriMCDplot',date],fighandle);



%% generate data for animation (gif) of MCD task
clear; clc

%%% set these in setParams.m
% params.armlengths(1) = 1.56*1.7 - 0.24;    % from kinematicChains.tex
% params.armlengths(2) = 1.3*1.7 - 0.24;     % " 
params = setParams;
Xmin = -4; Xmax = 5; Ymin = 1; Ymax = 5.5;

%%% Npts = 150; 
Npts = 500; 
locs = [(Xmax-Xmin)*rand(Npts,1)+Xmin,(Ymax-Ymin)*rand(Npts,1)+Ymin];
locs(:,1) = 2*locs(:,1) - (Xmin+Xmax)/2;
locs(:,2) = 2*locs(:,2) - (Ymin+Ymax)/2;
%%% for i = 1:Npts, fprintf('(%f,%f),\n',locs(i,1),locs(i,2)); end

% generate (right) shoulder and elbow angles
clear rRSA rREA
T1 = 20; T2 = 10; T3 = 5;
T = T1+T2+T3;
rRSA0 = 35;                         
rRSE0 = 70;
rRSA(1:T1) = rRSA0:1:(rRSA0+T1-1);  
rREA(1:T1) = rRSE0:2:(rRSE0+2*(T1-1));
rRSA((T1+1):(T1+T2)) = rRSA(T1):1:(rRSA(T1)+1*(T2-1));
rREA((T1+1):(T1+T2)) = rREA(T1):-2:(rREA(T1)-2*(T2-1));
rRSA((T1+T2+1):(T1+T2+T3)) = rRSA(T1+T2):-1:(rRSA(T1+T2)-1*(T3-1));
rREA((T1+T2+1):(T1+T2+T3)) = rREA(T1+T2):-1:(rREA(T1+T2)-1*(T3-1));

th = [rRSA',rREA']/180*pi + [-pi/2,0];
pos = FK2link(th,params.roboparams,0);
figure(13123);
plot(pos(:,1),pos(:,2));


% now compute "error vector" from hand to target
targetAngles = [60,85];
targetPos = FK2link(targetAngles/180*pi+[-pi/2,0],params.roboparams,0);
posError = pos - targetPos;
figure(13123); 
alp = 0.2;
locOffset = [0,0];
for t = 1:T
    clf; hold on;
    plot(pos(:,1),pos(:,2));
    scatter(locs(:,1)+locOffset(1),locs(:,2)+locOffset(2),'k.');
    axis([Xmin,Xmax,Ymin,Ymax]);
    scatter(pos(t,1),pos(t,2),'kx');
    scatter(targetPos(1),targetPos(2),'x');
    drawnow;
    pause(0.1);
    hold off;
    
    locOffset = locOffset - alp*posError(t,:);
    % locs = locs - alp*posError(t,:);
    % locs(:,1) = mod(locs(:,1)-Xmin,Xmax-Xmin)+Xmin;
    % locs(:,2) = mod(locs(:,2)-Ymin,Ymax-Ymin)+Ymin;
     
end


%% get movie of BMI reaches
clear; clc;

% load([getdir('data'),'RBMish/BMI/','wts_rEFH_spikecounts_160407_M1S1_CD1.mat'])
load([getdir('data'),'RBMish/BMI/','wts_rEFH_spikecounts_160627_M1_CD1.mat'])
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[R,X,Q] = params.getTestData(dataclass);
[~,Z] = updownRDBN(R,wts,params,Q.T);


%% see prev.
TOWRITE = 0;

% init
filename = 'point2pointReaching.gif';
past = 10;
k0 = 500 + past; % 1 + past;
kf = 200 + k0 - past; % size(trainData.Z,3); %%% probably want to shorten this
Xmin = min(X);
Xmax = max(X);
rmin = 0;                           
rmax = max(R);
H = figure(1);
colormap parula
Fs = params.Fs;
dt = params.Nmsperbin/1000;
setColors;

% for the circle plot
Nsegments = 32;
segs = 1:1:Nsegments;
th = 2*pi.*segs/(Nsegments) - pi;
xCirc = 0:0.01:1;
xslice = cos(th).*xCirc';
yslice = sin(th).*xCirc';
colorslices = @(iA,iB)(fill([xslice(:,iA);flipud(xslice(:,iB))],...
    [yslice(:,iA);flipud(yslice(:,iB))],1:11));
figure(3);
colormap([ones(11,1)  linspace(1,0,11)'  linspace(1,0,11)']);
iX = 1; iY = 2; iXdot = 3; iYdot = 4; iXddot = 5; iYddot = 6;


% loop through time
for k = k0:kf
    x = X((k-past):k,iX);
    y = X((k-past):k,iY);
    xdot = X((k-past):k,iXdot);
    ydot = X((k-past):k,iYdot);
    [phi,v] = cart2pol(X((k-past):k,iXdot),X((k-past):k,iYdot));
    r = R(k,:);
    z = Z(k,:);
    
    figure(1);
    
    subplot(2,2,1);
    squeezeplot(x,y);
    xlabel('x'); ylabel('y'); title('position')
    axis equal
    axis([Xmin(iX),Xmax(iX),Xmin(iY),Xmax(iY)]);
    
    subplot(2,2,2);
    squeezeplot(xdot,ydot);
    xlabel('$\dot x$','Interpreter','Latex'); 
    ylabel('$\dot y$','Interpreter','Latex'); 
    title('velocity')
    axis equal
    axis([Xmin(iXdot),Xmax(iXdot),Xmin(iYdot),Xmax(iYdot)]);
    
    subplot(2,2,3);
    %%%bar(1:length(r),r,1,'FaceColor',XTRAcolor,'EdgeColor','none');
    %%%xlabel('cell number'); ylabel('spikes/bin');
    %%%axis([1,length(r),rmin,rmax]);
    NN = floor(sqrt(length(r))); imagesc(reshape(r(1:NN^2),[NN,NN]));
    set(gca, 'XTick', []); set(gca, 'YTick', []); 
    title('neural responses');
    axis equal tight
    
    
    subplot(2,2,4);
    %%%bar(1:length(z),z,1,'FaceColor',EFHcolor,'EdgeColor','none');
    %%%axis([1,length(z),0,1]);
    %%%xlabel('unit number'); ylabel('p(spikes)'); title('rEFH responses');
    NN = floor(sqrt(length(z))); imagesc(reshape(z(1:NN^2),[NN,NN]));
    set(gca, 'XTick', []); set(gca, 'YTick', []); 
    title('rEFH responses');
    axis equal tight
   
    
    figure(3); clf; hold on;
    plot(xCirc,sqrt(1 - xCirc.^2),'k','linewidth',1.0);
    plot(xCirc,-sqrt(1 - xCirc.^2),'k','linewidth',1.0);
    plot(-xCirc,sqrt(1 - xCirc.^2),'k','linewidth',1.0);
    plot(-xCirc,-sqrt(1 - xCirc.^2),'k','linewidth',1.0);
    % plot(xslice,yslice,'k','linewidth',1.0);
    baseinds = diff(phi < th,[],2)*(1:(Nsegments-1))';
    indsA = mod(baseinds-1,32) + 1;
    indsB = mod(baseinds,32) + 1;
    colorslices(indsA,indsB);
    hold off;
    axis equal tight
   % xlabel('$\dot x$','Interpreter','Latex'); 
   % ylabel('$\dot y$','Interpreter','Latex'); 
   % title('velocity (polar coordinates)')
   
    
    if TOWRITE
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if k == k0;
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.05);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
        end
        
    else
        pause(dt);
    end
    
end





















































































