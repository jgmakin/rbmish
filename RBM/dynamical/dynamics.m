%% 
clear; clc; close all;

params = setParams;


%% for testing/debugging
tic
ldsDATA = getLDSdata(1000,params);
toc

%% plot the trajectories
[Ncases, Ndims, Nmods, T] = size(ldsDATA.S);

samples = 1:T;
% Ndims = params.Ndims;
NSind = strcmp(params.mods,params.NS);
smin = params.smin;
smax = params.smax;
srange = [params.smin(:,NSind) params.smax(:,NSind)];
thmin = params.thmin;
thmax = params.thmax;
th0 = get2DOutline(thmin,thmax,42); 
eTotal = zeros(Ncases,Ndims,T);


for iTraj = 1:Ncases
    
    thisS = squeeze(ldsDATA.S(iTraj,:,NSind,samples));
    thisX = params.dynamics.C*squeeze(ldsDATA.Z(iTraj,:,samples));
    thisShat = squeeze(ldsDATA.Y(iTraj,:,1,:));
    
    
    figure(45); clf;
    hold on;
    plot(thisS(1,:),thisS(2,:),'m')
    plot(thisX(1,:),thisX(2,:),'k')
    plot(th0(:,1),th0(:,2),'k')
    axis equal
    hold off;
    
    
    if 0
        for t = 1:T
            theseNoisyHills = displayshape(ldsDATA.R(iTraj,:,t),params);
            % figure(46); clf;
            % imagesc(theseNoisyHills{NSind});
            % axis image off; axis xy;
            figure(46), hold on;
            scatter(thisX(1,t),thisX(2,t),'b'); 
            scatter(thisShat(1,t),thisShat(2,t),'r'); 
            pause(0.05)
        end
    end
    
    
    % scatter them
    figure(47); clf;
    e = thisShat - thisX;
    scatter(e(1,:),e(2,:))
    
    % plot the correct Shat
    figure(45);
    hold on;
    plot(thisShat(1,:),thisShat(2,:),'c')
    hold off;
    
    mean(e');
    pause()
    
    eTotal(iTraj,:,:) = e;
    
end
eTotal = longdata(eTotal);
cov(eTotal)




%%

foo = squeeze(ldsDATA.R(iTraj,:,:));
for i=1:length(barr)
    figure(45); 
    scatter(thisShat(1,barr(i)),thisShat(2,barr(i)),'rx')
    
    caca = displayshape(foo(:,barr(i)),params);
    %%%% replace me
    wrappeddecode(caca{1},params);
    %%%%
end



%%
% make sure things work for uniform stimuli, at least
clc
S = getUniformStims(15,params);
params.Ncases = 15;
[D0,S,tilde,Q] = DATAGENPP(size(S,1)/params.Ncases,params,'stimuli',S);




figure(101); clf; hold on;
Shat = zeros(size(S));
for i = 1:size(D0,1)
    for j = 1:size(D0,3)
        TT = displayshape(D0(i,(params.t+1):end,j),params);
        %         figure(46); clf;
        %         imagesc(TT{1});
        %         axis image off; axis xy;
        %         pause(0.05)
        Shat(i,:,1,j) = decode(TT{1},[params.smin(:,1) params.smax(:,1)],params,'torusCoM')';

    end
end
figure(101); clf; hold on;
foo = longdata(S);
scatter(foo(:,1),foo(:,2));
foo = longdata(Shat);
scatter(foo(:,1),foo(:,2),'m');
E = Shat - S;
foo = longdata(E);
caca = longdata(S);
hh = figure(102); clf;
quiver(caca(:,1),caca(:,2),foo(:,1),foo(:,2));
max(abs(E(:)))



%% get the relative covariances!!
popind = 2;                         % the PPC is in prop coordinates
A = params.dynamics.A;
G = params.dynamics.G;
C = params.dynamics.C;
SigmaX = params.dynamics.SigmaX;
tuningCov = computetuningcovs(params);
R = longdata(D0(:,(params.N^2+1):end,:));
eta = squeeze(sum(R,2));
invTCprop = inv(tuningCov{popind});
ExpCovEgivenR = inv(invTCprop)*mean(1./eta);

[Pthe,L,G] = dare(A',C',G*SigmaX*G',ExpCovEgivenR);

inv(C'*inv(ExpCovEgivenR)*C + inv(Pthe))
stats{2}{2}.cov
% Pemp = (A*inv(C'*inv(stats{2}{2}.cov)*C)*A' + B*SigmaX*B');


%% for testing/debugging
tic
N = params.N;
setColors;


% give the KF all the parameters it needs
KFparams.T = T;
KFparams.A = params.dynamics.A;
KFparams.C = params.dynamics.C;
KFparams.SigmaX = params.dynamics.SigmaX;
KFparams.mu0 = [params.dynamics.muX0; params.dynamics.muV0];
KFparams.Info0 = setInfoMatrix(params.dynamics.SigmaX0,...
    params.dynamics.SigmaV0,params.Ndims);


% run the Kalman filter
%%%%%%%%%
[pKF, pSENSORY] = KF4PPC(D0(:,(end-N^2+1):end,:),KFparams,'some name');
%%%%%%%%% this is now the wrong syntax


pSENSORY.color = PROPcolor;     pSENSORY.name = 'raw';
pKF.color = OPTcolor;           pKF.name = 'opt';

toc

%% 
t0 = 50;
j = ceil(rand*params.Ncases);

stats = dispFilterErrCovs(t0,S0,params,pSENSORY,pKF);
plotFilterErr(t0,1:size(S0,1),S0,params,pSENSORY,pKF);


























