function [Zout,Uout,LDSparams] = getCenterOutReacher(T,params)
%%%%%%%%%%% you don't really need to return X, probably
%
% Load in some of H.H.Shen's center-out reach data (Dmitri), fit a LDS,
% then....

%-------------------------------------------------------------------------%
% Revised: 03/13/15
%   -functionized, rationalized, etc.
% Created: 03/12/15
%   -by JGM
%-------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%
%%%% TO DO
% (4) printf warning if psoitions or velocities leave the feasible space??
% (7) ok, maybe not clip all trajectories to the same size, but you should
% at least get rid of the ones that don't make it anywhere near the target.
%%%%%%%%%%%%%%%%%%


if strcmp(params.machine,'domestica')
    yrclass = 'gpuArray';
else
    yrclass = 'double';
end


% useful params
HHSdt = 1/240;                              % fact
NsamplesPerBin = params.dynamics.dt/HHSdt;  % convert to longer bins

% Ns
Ncases = params.Ncases;
Ndims = params.Ndims;
LDSparams.Ndims = Ndims;
LDSparams.Nstates = Ndims*2;                 % 2nd-order
LDSparams.Noutputs = Ndims;                  % control position(s)
LDSparams.Ninputs = LDSparams.Noutputs;       % to stand a chance at control!
LDSparams.cntr = params.dynamics.muX0;


% do things!!
[S,targets] = getHHSdata('D080620');
[X,endinds,LDSparams] = fitStateTransitionModel(S,HHSdt,NsamplesPerBin,...
    LDSparams,yrclass);
Z = arrayfun(@(iTraj)(lengthenTrajectories(X,endinds,targets,0,T)),...
    1:Ncases,'UniformOutput',false);
Z = shiftdim(cat(3,Z{:}),2);
[LDSparams.muX,LDSparams.SigmaX] = scaleDownTransisionNoise(Ndims,...
    LDSparams.muX,LDSparams.SigmaX);
[Zout,Uout,LDSparams.G] = getCenterOutControls(Z,LDSparams,...
    params.smin,params.smax,params.mods);

% plot autocorrelations
% fignum = plotAutocorrelations(Z(:,1:Ndims,:),...
%     smin(:,1),smax(:,1),LDSparams.dt,fignum);
% fignum = plotAutocorrelations(Z(:,Ndims+(1:Ndims),:),...
%     smin(:,2),smax(:,2),LDSparams.dt,fignum);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [St,targetLocsT] = getHHSdata(tag)

datadir = '../HHS/extracteddata/';
load([datadir,'KFtuningdataHHS',tag]); % "the tuning series"
%%% NB: this assumes that this data set has been created.  If it hasn't,
%%% use extractKFdata.m (in the HHS directory). NOTE THAT THIS M-FILE WILL
%%% SAVE THE .MATs IN YOUR CURRENT DIR, *NOT* \extracteddata\.
fprintf('\n\nLoaded data...\n\n');
%%% You could also try JEO's data:
% load('../../#DATA/Indy_datafiles/HandControl_indy_bmi5_20140505_1.mat')
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [X,endinds,LDSparams] = fitStateTransitionModel(...
    S,dt,NsamplesPerBin,LDSparams,yrclass)
% fit A, SigmaX, muX

% Ns
Ntrials = length(S);
Nstates = LDSparams.Nstates;
Ndims = LDSparams.Ndims;
binsize = NsamplesPerBin*dt;

% store the trial-ending indices (but starting with 0)
ind = 0;
endinds = zeros(Ntrials+1,1,yrclass);
for iTrial = 1:Ntrials
    edges = S(iTrial).t(1):binsize:S(iTrial).t(end);
    Nbins = length(edges) - 1;
    ind = ind+Nbins;
    endinds(iTrial+1) = ind;
end

% now store the data (you make two separate loops for malloc purposes)
X = zeros(ind,Nstates,yrclass);
for iTrial = 1:Ntrials
    switch Nstates/Ndims %%% not very elegant
        case 1
            thisX = [S(iTrial).pos];
        case 2
            thisX = [S(iTrial).pos S(iTrial).vel];
        case 3
            thisX = [S(iTrial).pos S(iTrial).vel S(iTrial).acc];
    end
    Nbins = length(S(iTrial).t(1):binsize:S(iTrial).t(end))-1;
    X((endinds(iTrial)+1):(endinds(iTrial+1)),:) = squeeze(mean(reshape(...
        thisX(1:(Nbins*NsamplesPerBin),:)',...
        [Nstates,NsamplesPerBin,Nbins]),2))';
end



if 1
    %%%%%%%%%% this assumes x and v are both present
    %%%% this should be cleaned up
    % you set these
    speedThr = 120;
    % speedThr = 60;
    distThr = 40; % cm
    
    cntr = LDSparams.cntr;
    x = X(:,1:Ndims);
    cntrDist = squeeze(sqrt(sum(bsxfun(@minus,x,shiftdim(cntr,-1)).^2,2)));
    v = X(:,(Ndims+1):2*Ndims);
    mvmtSpeed = sqrt(sum(v.^2,2));
    inds = gatherCenterOutInds(cntrDist,distThr,mvmtSpeed,speedThr,...
        size(X,1),0,[],[],[]);
    
    allinds = arrayfun(@(iReach)(inds(1,iReach):inds(2,iReach)),...
        1:size(inds,2),'UniformOutput',false);
    futureInds = arrayfun(@(iReach)((inds(1,iReach)+1):inds(2,iReach)),...
        1:size(inds,2),'UniformOutput',false);
    pastInds = arrayfun(@(iReach)(inds(1,iReach):(inds(2,iReach)-1)),...
        1:size(inds,2),'UniformOutput',false);
    
    Xf = X(cat(2,futureInds{:}),:);
    Xp = X(cat(2,pastInds{:}),:);
    X = X(cat(2,allinds{:}),:);
    
    endinds = [0,cumsum(arrayfun(@(iReach)(length(allinds{iReach})),...
        1:length(allinds)))];
    fprintf('\nFraction removed: %.3f\n\n',(size(x,1)-size(X,1))/size(x,1));
    %%%%%%%%%%

else
    Xf = X(~ismember(1:size(X,1),endinds+1),:);
    Xp = X(~ismember(1:size(X,1),endinds),:);
end


% regress
[beta,RsqCV,ResCV] = linrgsLOOCV(Xp,Xf);
fprintf('A fit with R^2 = %0.3f\n\n',RsqCV);

% store
LDSparams.A = beta';
LDSparams.muX = mean(ResCV)';
LDSparams.SigmaX = cov(ResCV);
LDSparams.dt = binsize;

% exit
fprintf('\n\nfit dynamics (A, SigmaX, muX)...\n\n');

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [muX,SigmaX] = scaleDownTransisionNoise(Ndims,muX,SigmaX)
% assume a lot of the state-transition "noise" can be explained by a ctrl

% integers
Nstates = size(SigmaX,2);
inds = (Nstates-Ndims+1):Nstates;

% scale
muX = muX/1000;
SigmaX = SigmaX/1000;
% SigmaX(inds,:) = SigmaX(inds,:)/10;
% SigmaX(:,inds) = SigmaX(:,inds)/10;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Zout,Uout,G] = getCenterOutControls(Z,LDSparams,smin,smax,mods)
% find some controls that make it do center-out reaches!


% Ns
Ninputs = LDSparams.Ninputs;
Noutputs = LDSparams.Noutputs;
Nstates = LDSparams.Nstates;
Ndims = LDSparams.Ndims;
dt = LDSparams.dt;

% x[t+1] = A*x[t] + G*u[t] + noise,     y[t] = C*x[t] + noise
m = 5;
G = zeros(Nstates,Ninputs);
G((end-Ninputs+1):end,(end-Ninputs+1):end) = eye(Ninputs)*dt/m;
%%% must match what's in params.m!!
C = zeros(Noutputs,Nstates);
C(1:Ndims,1:Ndims) = eye(Ndims);
LDSparams.G = G;
LDSparams.C = C;

% track the "true" monkey reaches
Ref = Z(:,1:Ndims,:);
umin = smin(:,strcmp(mods,'Efference-Copy'));
umax = smax(:,strcmp(mods,'Efference-Copy'));
[Zout,Uout] = generateTrackingControls(Ref,umin,umax,LDSparams);




end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function longTraj = lengthenTrajectories(X,endinds,targets,TOPLOT,T)
% Turn center-out trajectories into center-out-center-out etc. etc.
%
% This could be sped up through vectorization....


% Ns
Ntargets = size(targets,1);
Nstates = size(X,2);

% things you need in the loop
cntr = mean(X(endinds(1:end-1)+1,:));
trajTargs = assignTrajsToTargets(X(endinds(2:end),:),targets,TOPLOT,4121);


% init
longTraj = zeros(Nstates,T,'like',X);
ind = 1;
NEEDMOREBATCHES = true;

% while there are fewer than T batches....
while NEEDMOREBATCHES

    % check how much room is left in longTraj
    remainingSamples = T-ind+1;
    
    % get a random target (so that an equal number of each targ is chosen)
    thisTarg = ceil(Ntargets*rand);
    trajsToThisTarg = find(trajTargs == thisTarg);
    thisTrial = trajsToThisTarg(ceil(length(trajsToThisTarg)*rand));
    thisCenterOut = X((endinds(thisTrial)+1):endinds(thisTrial+1),:);
    
    
    % now simulate a reach back to the center
    R = [cos(pi), -sin(pi); sin(pi), cos(pi)];
    R = blkdiag(R,R);
    xF = thisCenterOut(end,:);
    thisOutCenter = bsxfun(@plus,(R*bsxfun(@minus,thisCenterOut,xF)')',cntr);
    
    
    % concatenate together
    %%%% hold briefly??
    inAndOut = cat(1,thisCenterOut,thisOutCenter);
    trajLength = size(inAndOut,1);
    
    
    % check if there's room for this traj
    if trajLength > remainingSamples
        NEEDMOREBATCHES = false;
        longTraj(:,ind:end) = inAndOut(1:remainingSamples,:)';
    else
        longTraj(:,ind:(ind+trajLength-1)) = inAndOut';
        ind = ind+trajLength;
    end
    
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Xclipped,Zclipped,Uclipped,fignum] =...
    clipTrajectories(Xout,Zout,Uout,endinds,fignum)
% make controls and trajectories of a single length!


% Ns
trialLengths = diff(endinds);           % in samples
figure(fignum);
hist(trialLengths);                     % ...looks like 18 is good number 
title('N samples per trial')            %   of samples to use
Msamples = 18;                          %   
fignum = fignum+1;

startindsLongTrials = endinds([(trialLengths > Msamples);false])+1;
Mtrials = length(startindsLongTrials);

Xclipped = arrayfun(@(ind)(Xout(startindsLongTrials+ind,:)),1:Msamples,...
    'UniformOutput',false);
Xclipped = cat(3,Xclipped{:});
Zclipped = arrayfun(@(ind)(Zout(startindsLongTrials+ind,:)),1:Msamples,...
    'UniformOutput',false);
Zclipped = cat(3,Zclipped{:});
Uclipped = arrayfun(@(ind)(Uout(startindsLongTrials+ind,:)),1:Msamples,...
    'UniformOutput',false);
Uclipped = cat(3,Uclipped{:});

%%%%%%%%%%
%%%% change to use Ndims etc. rather than 1,2,3,4?
%%%%%%%%%%
figure(fignum); clf; figure(fignum+1); clf;
for iTrial = 1:Mtrials
    
    % plot true trajs and the ones you get via the controls
    figure(fignum);
    hold on;
    plot(squeeze(Xclipped(iTrial,1,:)),squeeze(Xclipped(iTrial,2,:)),'k')
    plot(squeeze(Zclipped(iTrial,1,:)),squeeze(Zclipped(iTrial,2,:)),'g')
    hold off;
    
    % ditto for velocities!
    figure(fignum+1);
    hold on;
    plot(squeeze(Xclipped(iTrial,3,:)),squeeze(Xclipped(iTrial,4,:)),'k')
    plot(squeeze(Zclipped(iTrial,3,:)),squeeze(Zclipped(iTrial,4,:)),'g')
    hold off;
end

% entitle the figures
figure(fignum);
titlestr = sprintf('Clipped trajectories (to %i samples)',Msamples);
title(titlestr);
axis equal

figure(fignum+1);
titlestr = sprintf('Clipped velocities (to %i samples)',Msamples);
title(titlestr);
axis equal

fignum = fignum+2;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function closestTarg = assignTrajsToTargets(Zfinal,targets,TOPLOT,clippedTrajFignum)

% Ns
Ntrials = size(Zfinal,1);
Ndims = size(targets,2);

% find the (nominal) targets closest to the ending point of each traj
vecsToTargs = bsxfun(@minus,Zfinal(:,1:Ndims),shiftdim(targets',-1));
distsToTargs = shortdata(Ntrials,3,...
    sum(longdata(vecsToTargs).*longdata(vecsToTargs),2));
[~,closestTarg] = min(distsToTargs,[],3);


% check
if TOPLOT
    figure(clippedTrajFignum); hold on;
    map = colormap;
    for iTarg = 1:size(targets,1)
        foo = Zfinal((closestTarg==iTarg),1:Ndims);
        scatter(foo(:,1),foo(:,2),[],map(iTarg*8,:),'o')
        scatter(targets(iTarg,1),targets(iTarg,2),[],map(iTarg*8,:),'x');
    end
    hold off;
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function fignum = plotAutocorrelations(X,xmin,xmax,dt,fignum)

% Ns
Ntrials = size(X,1);
tauMax = 12;

% init
acAngle = 0;      
acCartesian = 0;    
acCos = 0;     
acSin = 0;
foo = 0;

% wrap the data
Xscaled = shortdata(Ntrials,3,scalefxn(squeeze(longdata(X)'),...
    xmin,xmax,zeros(2,1),2*pi*ones(2,1))');

% loop through trials (trajectories)
for iTrial = 1:Ntrials
    
    % autocorrelation of the angle
    x = squeeze(X(iTrial,1,:));
    y = squeeze(X(iTrial,2,:));
    ang = atan2(y,x);
    acAngle = acAngle + xcorr(ang,ang,tauMax,'coeff');
    
    % autocorrelation of the variable itself
    acCartesian = acCartesian + xcorr([x y],tauMax,'coeff'); 
    
    % circular autocorrelation of wrapped version
    xscaled = squeeze(Xscaled(iTrial,:,:))';
    acCos = acCos + xcorr(cos(xscaled),tauMax,'none');
    acSin = acSin + xcorr(sin(xscaled),tauMax,'none');
    
    % foo!
    caca = atan2(xscaled(:,2),xscaled(:,1));
    foo = foo + xcorr(cos(caca),cos(caca),tauMax,'none');
    
    
end
acAngle = acAngle/Ntrials;
acCartesian = acCartesian/Ntrials;
acCos = acCos/Ntrials;
acSin = acSin/Ntrials;

% plot
figure(fignum); clf;
hold on;
plot((-tauMax:tauMax)*dt,acAngle);
plot((-tauMax:tauMax)*dt,acCartesian(:,1));
plot((-tauMax:tauMax)*dt,acCartesian(:,4));
plot((-tauMax:tauMax)*dt,acCos(:,1)/max(acCos(:,1))); %%%
plot((-tauMax:tauMax)*dt,acCos(:,4)/max(acCos(:,4))); %%%
plot((-tauMax:tauMax)*dt,acSin(:,1)/max(acSin(:,1))); %%%
plot((-tauMax:tauMax)*dt,acSin(:,4)/max(acSin(:,4))); %%%
plot((-tauMax:tauMax)*dt,foo/max(foo));
legend('angle','Cartesian x-coord','Cartesian y-coord',...
    'cos(x)','cos(y)','sin(x)','sin(y)','foo!');
title('Autocorrelations');
hold off;
fignum = fignum + 1;

% exit
fprintf('\n\ncomputed autocorrelations...\n\n');


end
%-------------------------------------------------------------------------%







%-------------------------------------------------------------------------%
function dasfdasfsdaf(Xout)
%Do it a different way: 

%%%%%%%%%
% What you really need/ought to do to complete this is modify all your
% Kalman filter files to work on sets of trajectories *of different
% (temporal) lengths*!!
%%%%%%%%%


% fit a six-state model directly, forcing the third-order stuff to be a
% noisily observed state!
X = Xout';
Xf = Xout(~ismember(1:size(Xout,1),endinds+1),:)';
Xp = Xout(~ismember(1:size(Xout,1),endinds),:)';
x00 = Xout(endinds(1:end-1)+1,:)';


Nsamples = size(X,2);
NN = Nsamples - Ntrials; % you remove either the first or the last

% sufficient stats: first-order
xpctT.mu0 = mean(x00,2);

% sufficient stats: second-order
xpctT.x0x0 = x00*x00'/Ntrials;                      % ...the first
xpctT.XpXp = Xp*Xp'/NN;                             % ...all but the last
xpctT.XfXp = Xf*Xp'/NN;                             % ...consecutive states
xpctT.XfXf = Xf*Xf'/NN;                             % ...all but the first
xpctT.XX = X*X'/Nsamples;                           % ...all states

Nx = 4;
Nu = 2;


% muU = (<U0> + (T-1)<Uf>)/T 
muU = (xpctT.mu0(Nx+(1:Nu),1) + (T-1)*xpctT.XfXp(Nx+(1:Nu),end))/T;
xpctT.mu0(Nx+(1:Nu),1) = muU;       
xpctT.XfXp(Nx+(1:Nu),end) = muU;

% UU = (<U0*U0> + (T-1)<Uf*Uf>)/T, 
UU = (xpctT.x0x0(Nx+(1:Nu),Nx+(1:Nu)) +...
    (T-1)*xpctT.XfXf(Nx+(1:Nu),Nx+(1:Nu)))/T;
xpctT.x0x0(Nx+(1:Nu),Nx+(1:Nu)) = UU; 
xpctT.XfXf(Nx+(1:Nu),Nx+(1:Nu)) = UU;

% zeros in Gamma
Uf = xpctT.XfXp(Nx+(1:Nu),end);
Zp = xpctT.XpXp(1:end-1,end);
xpctT.XfXp(Nx+(1:Nu),1:end-1) = Uf*Zp';
%%% xpctT.XfXp/xpctT.XpXp

% zeros in SigmaZ
XfbarUbar = xpctT.XfXp(1:Nx,end)*xpctT.XfXp(Nx+(1:Nu),end)';
xpctT.XfXf(1:Nx,Nx+(1:Nu)) = XfbarUbar;
xpctT.XfXf(Nx+(1:Nu),1:Nx) = XfbarUbar';
%%% xpctT.XfXf - xpctT.XfXp/xpctT.XpXp*xpctT.XfXp'

end
%-------------------------------------------------------------------------%

