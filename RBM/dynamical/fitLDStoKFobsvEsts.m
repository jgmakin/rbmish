function LDSparams = fitLDStoKFobsvEsts(LDSdataTest,model,params,varargin)
% fitLDStoKFobsvEsts    Fit a LDS using the outputs of a secret KF
%
% It may happen (strangely) that you have access to a Kalman filter's
% prediction of the emission, yhat, but not the parameters of the linear
% dynamical system assumed by the KF---in your case, because the "KF" is
% really a recurrent, exponential-family harmonium.
% 
% This function takes as input LDSdata (see getLDSdata.m), and the weights
% and parameters of a model: either an rEFH (wts,params) or a linear
% dynamical system used for Kalman filtering (LDSparams,params)---probably
% one learned with EM.  It initializes the LDSparams to be learned,
% {A,C,SigmaZ}, according to the last argument; and then adjusts them by
% descending the gradient of the squared error between the model's decoded
% outputs and the yhats produced by a KF that uses the current parameters.
% This gradient must be computed with a backwards recursion on the KF
% recursions, which is done inside BPTT4KF.m.
%
% The final argument controls the initialization of the parameters:
%   'true'      -the *true* params of the underlying LDS
%   'random'    -random parameters
%   'EM2'       -parameters learned with EM (2nd-order)
%   LDSparams   -these LDS parameters
%
%   USAGE:
%{
        load('dynamical\finalwts\wts1DrEFHManyXprmts.mat')
        LDSdata = getLDSdata(params);
        LDSparams = fitLDStoKFobsvEsts(LDSdata,Allwts{1},params,'true');

        load dynamical\finalwts\LDSparamsEM2ndOrd1DrEFHManyXprmts.mat
        LDSdata = getLDSdata(params);
        LDSparams = fitLDStoKFobsvEsts(LDSdata,Allparams(1),params,'true');
%}

%-------------------------------------------------------------------------%
% Revised: 06/26/15
%   -added arguments for learning rate and initialization
% Revised: 06/25/15
%   -moved training-data error calculations into BPTT4KF
%   -modified to allow for EM-based models, as well as rEFHs
% Created: 06/24/15
%   by JGM
%-------------------------------------------------------------------------%

%%%%%
% Currently broken because of the new KF4PPC and LDSdata.... (01/08/17)
%%%%%

% init
Nepochs = 50;
plotHandle = initPlots(2014);
fignum = 4;
figure(fignum); clf;
EMBASED = isstruct(model);
T = 1000;
Ntraj = 40;

% targets for *testing* data (fixed)
YSTARtest = getKFtargets(model,LDSdataTest,EMBASED,params);
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

% malloc
MSE = zeros(Nepochs,1);

% train for Nepochs
for iEpoch = 1:Nepochs
    
    % make new training data every epoch
    LDSdataTrain = getLDSdata(Ntraj,T,dataclass,params);
    if iEpoch==1
        if isstruct(varargin{1})
            LDSparams = varargin{1};
            mass0 = 96095; % 5000; % 200;
        else
            [LDSparams,mass0] = initLDSparams(class(LDSdataTrain.Y),...
                varargin{1},params);
        end
    end
    
    % filter it with the MODEL; these will be the training targets
    YSTARtrain = getKFtargets(model,LDSdataTrain,EMBASED,params);
    
    % get the errors of the model *on testing data*
    MSE(iEpoch) = getKFmirrorErr(LDSdataTest,YSTARtest,LDSparams,0);
    
    % print errors
    fprintf('epoch: %i; terror: %2.4d\n',iEpoch,MSE(iEpoch));
    set(0,'CurrentFigure',figure(2014));
    hold on;
    set(plotHandle(1),'XData',1:iEpoch,'YData',MSE(1:iEpoch),'color','r');
    hold off;
    drawnow
    
    % use posterior means from the rEFH to learn LDS parameters, w/BPTT
    mAC = 1.03^iEpoch*mass0; % params.massUpdate(mass0,iEpoch);
    LDSparams = BPTT4KF(YSTARtrain,LDSdataTrain,LDSparams,mAC);
    
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function plotHandle = initPlots(fignum)

figure(fignum); % clf;
hold on;
plotHandle(1) = plot(NaN,NaN);
hold off;

% subplot(1,2,1); hold on; 
% plotHandle(1) = plot(NaN,NaN);
% hold off;
% subplot(1,2,2); hold on;
% plotHandle(2) = plot(NaN,NaN);
% hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [LDSparams,mass0] = initLDSparams(yrclass,INITTYPE,params)

fprintf('initializing parameters...\n');
switch INITTYPE
    case 'random'
        [Ny,Nx] = size(params.dynamics.C);
        LDSparams.A = randn(Nx,Nx,yrclass);
        %     if isfield(LDSdataTrain,'U'),
        %         LDSparams.B = randn(Nx,Nu,yrclass);
        %     end
        foo = randn(Nx,Nx,yrclass)/500;
        %%%
        LDSparams.SigmaX = params.dynamics.SigmaX; %%%foo'*foo;
        %%%
        LDSparams.muX = zeros(Nx,1); %%% randn(Nx,1,'like',LDSparams.Y);
        LDSparams.C = randn(Ny,Nx,yrclass);
        %%%
        foo = randn(Nx,Nx,yrclass);
        LDSparams.Info0 = [0 0; 0 2e9]; %%% foo'*foo;
        LDSparams.mu0 = zeros(Nx,1,yrclass); %%% randn(Nx,1,yrclass);
        %%%
        mass0 = 100;
    case 'true'     % LDS params of the system originally learned
        LDSparams = getLDSparams(params.dynamics);
        mass0 = 250;
    case 'EM2'      % LDS params learned by EM, second-order
        load([getdir('data'),'RBMish/EMparams/LDSOrd2_1D_LTI-PPC_ManyXprmts.mat'])
        LDSparams = Allparams(2);
        mass0 = 500;
    otherwise
        error('whoops')
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function MSE = getKFmirrorErr(LDSdata,YSTAR,LDSparams,fignum)

LDSparams.SigmaYX = LDSdata.SigmaYX;
pMirror = KF4PPC(LDSdata,LDSparams,'mirror');
MSE = mean((pMirror.Xpct(:) - YSTAR(:)).^2);

inds = 401:460;
figure(1); clf; hold on;
plot(squeeze(LDSdata.S(1,:,:,inds)))
plot(squeeze(pMirror.Xpct(1,inds)))
plot(squeeze(YSTAR(1,:,inds)))
hold off;
drawnow


if fignum
    figure(fignum); hold on; drawnow;
    plot(mean((pMirror.Xpct - YSTAR).^2,3));
    hold off
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function YSTAR = getKFtargets(model,LDSdata,EMBASED,params)

if EMBASED
    pMODEL = KF4PPC(LDSdata,model,'EM');
else
    [~,~,pMODEL] = EFHfilter(LDSdata,model,params);
end

YSTAR = pMODEL.Xpct;

end
%-------------------------------------------------------------------------%

