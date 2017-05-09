function [Rsqs,NdataTest] = filtersForNeuralData(wts,params,LDSparamsEM,...
    decodernames,TOPLOT)
% filtersForNeuralData
% USAGE:
%   load([getdir('data'),'RBMish/BMI/wts_rEFH_spikecounts_160407_M1S1_CD1.mat'],'wts','params');
%   [Rsqs,NdataTest] = filtersForNeuralData(wts,params,LDSparamsEM,decodernames,TOPLOT);
%
% Fit filters (and other decoders) to neural data

%-------------------------------------------------------------------------%
% Revised: 04/15/17
%   -made the various decoders requestable
% Created: 04/??/16
%   by JGM
%-------------------------------------------------------------------------%


% init
if checkGPUavailability, dataclass = 'gpuArray'; else, dataclass = 'double'; end
tplotmax = 60; % s
if TOPLOT, tvec = 1:ceil(tplotmax/(params.Nmsperbin/1000)); else, tvec=[]; end

% get data
[Rtrain,Xtrain,Qtrain] = getSingleTrainingSequence(dataclass,params);
[Rtest,Xtest,Qtest] = params.getTestData(dataclass);
SStot = sum((Xtest - mean(Xtest)).^2);
NdataTest = size(Rtest,1);
if ~isempty(tvec), PosVelPlot(tvec,Xtest,1); end




% (1) (static) linear regression
if any(strcmp(decodernames,'static'))
    [beta,Xhat,RsqStatic] = staticLinearDecoder(Rtrain,Xtrain,...
        Rtest,Xtest,SStot,tvec);
    Rsqs(:,:,:,strcmp(decodernames,'static'))       = RsqStatic;
end

% (2) observed LDS + KF
if any(strcmp(decodernames,'kfobs'))
    [LDSparamsObs, KFdstrbsObs, RTSSdstrbsObs, RsqKF_Obs, RsqRTSS_Obs] = ...
        fullyObservedLDSDecoder(Rtrain,Xtrain,Qtrain.T,Rtest,Xtest,SStot,tvec);
    Rsqs(:,:,:,strcmp(decodernames,'kfobs'))        =...
        RsqKF_Obs(:,:,:,any(strcmp(decodernames,'kfobs')));
    Rsqs(:,:,:,strcmp(decodernames,'rtssobs'))      =...
        RsqRTSS_Obs(:,:,:,any(strcmp(decodernames,'rtssobs')));
end

% (3,4) KF with latent state, plus regresion to kinematics
if any(ismember(decodernames,{'kfemstatic','kfemdynamic'}))
    [LDSparamsEM, Bxz_EM, LDSparamsObs, RsqKF_EMstatic,RsqKF_EMdynamic] =....
        latentStateLDSDecoder(Rtrain,Xtrain,Rtest,Xtest,params,LDSparamsEM,...
        SStot,tvec);
    Rsqs(:,:,:,strcmp(decodernames,'kfemstatic'))   =...
        RsqKF_EMstatic(:,:,:,any(strcmp(decodernames,'kfemstatic')));
    Rsqs(:,:,:,strcmp(decodernames,'kfemdynamic'))  =...
        RsqKF_EMdynamic(:,:,:,any(strcmp(decodernames,'kfemdynamic')));
end

% (5) UKF
if any(strcmp(decodernames,'ukf'))
    %if any(strncmp(decodernames,'ukf',3))
    %ukfname = decodernames{strncmp(decodernames,'ukf',3)};
    %Ntapspast = str2double(ukfname(4));
    %Ntapsfuture = str2double(ukfname(5));
    [UKFparams, RsqUKF] = UKFDecoder(Rtrain,Xtrain,Qtrain.T,...
        Rtest,Xtest,SStot,tvec);
    %%%Rsqs(:,:,:,strncmp(decodernames,'ukf',3))       = RsqUKF;
    Rsqs(:,:,:,strcmp(decodernames,'ukf'))          = RsqUKF;
end

% (6,7) rEFH
if any(ismember(decodernames,{'refhstatic','refhdynamic',...
        'refhstatic_stdnrml','refhdynamic_stdnrml'}))
    [RsqREFHstatic,RsqREFHdynamic] = linearREFHDecoder(Rtrain,Xtrain,Qtrain,...
        Rtest,Xtest,Qtest,wts,params,SStot,tvec);
    Rsqs(:,:,:,strcmp(decodernames,'refhstatic'))   =... 
        RsqREFHstatic(:,:,:,any(strcmp(decodernames,'refhstatic')));
    Rsqs(:,:,:,strcmp(decodernames,'refhdynamic'))  =... 
        RsqREFHdynamic(:,:,:,any(strcmp(decodernames,'refhdynamic')));
    Rsqs(:,:,:,strcmp(decodernames,'refhstatic_stdnrml'))   =... 
        RsqREFHstatic(:,:,:,any(strcmp(decodernames,'refhstatic_stdnrml')));
    Rsqs(:,:,:,strcmp(decodernames,'refhdynamic_stdnrml'))  =... 
        RsqREFHdynamic(:,:,:,any(strcmp(decodernames,'refhdynamic_stdnrml')));
end



% (8,9) Poisson emissions (particle filter)
if any(ismember(decodernames,{'particlefilter','particlesmoother'}))
    [LDSparamsPE, PFdstrbs, PSdstrbs, RsqPF, RsqPS] = ...
        fullyObservedPoissonDecoder(trainData,testData,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'particlefilter'))   =... 
        RsqPF(:,:,:,any(strcmp(decodernames,'particlefilter')));
    Rsqs(:,:,:,strcmp(decodernames,'particlesmoother')) =... 
        RsqPS(:,:,:,any(strcmp(decodernames,'particlesmoother')));
end

% (10,11) observations linear in polar transform of vel (particle filter)
if any(ismember(decodernames,{'polarpf','polarps'}))
    [LDSparamsPolar, PolarPFdstrbs, PolarPSdstrbs, RsqPolarPF, RsqPolarPS] = ...
        fullyObservedPolarDecoder(LDSparamsObs,trainData,testData,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'polarpf'))      =... 
        RsqPolarPF(:,:,:,any(strcmp(decodernames,'polarpf')));
    Rsqs(:,:,:,strcmp(decodernames,'polarps'))      =...
        RsqPolarPS(:,:,:,any(strcmp(decodernames,'polarps')));
end





if 0
    % Next-frame predictions
    YHATOBS = LDSparamsObs.C*KFdstrbsObs.XHATTU + LDSparamsObs.muYX;
    [~,Z0] = updownRDBN(Rtest,wts,params,Qtest.T);
    Z0 = shortdata(floor(size(Z0,1)/Qtest.T),3,Z0);
    tic
    params.Ncases = 1; %%% what is this supposed to do?
    YHATREFH = longdata(forwardGenerate(Z0,1,0,'inferredhidmeans',wts,params))';
    params.Ncases = 4;
    toc
    YHATEM = LDSparamsEM.C*KFdstrbsEM.XHATTU + LDSparamsEM.muYX;
    SStotY = sum((Rtest(2:end,:) - mean(Rtest(2:end,:))).^2);
    
    
    keyboard
    RsqsAvg = averageDecentRsqs(sum(diff(Rtest).^2),SStotY,3241)
    RsqsAvg = averageDecentRsqs(sum((YHATOBS(:,2:end)' - Rtest(2:end,:)).^2),SStotY,3242)
    RsqsAvg = averageDecentRsqs(sum((YHATEM(:,2:end)' - Rtest(2:end,:)).^2),SStotY,3244)
    RsqsAvg = averageDecentRsqs(sum((YHATREFH(:,2:end)' - Rtest(2:end,:)).^2),SStotY,3243)
end




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function PosVelPlot(tvec,X,CLEARPLOTS)

if usejava('desktop')
    % plot
    figure(101);
    if CLEARPLOTS, clf; end
    hold on
    plot(tvec,X(tvec,2))
    hold off;
    
    figure(102);
    if CLEARPLOTS, clf; end
    hold on
    plot(tvec,X(tvec,5))
    hold off;
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [beta, Xhat, RsqStatic] = staticLinearDecoder(Ytrain,Xtrain,...
    Ytest,Xtest,SStot,tvec)

% decoder
beta = linregress([Ytrain,ones(size(Ytrain,1),1)],Xtrain);
Xhat = [Ytest, ones(size(Ytest,1),1)]*beta;
SSerr = sum((Xhat - Xtest).^2);
RsqStatic = 1 - SSerr./SStot;

% plot?
if ~isempty(tvec), PosVelPlot(tvec,Xhat,0); end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [LDSparamsOBS, KFdstrbs, RTSSdstrbs, RsqKF, RsqRTSS] = ...
    fullyObservedLDSDecoder(Ytrain,Xtrain,T,Ytest,Xtest,SStot,tvec)

% Ns
Nexamples = size(Ytrain,1);
Ntraj = floor(Nexamples/T);

% train
LDSparamsOBS = learnfullyobservedLDS(shortdata(Ntraj,3,Ytrain),...
    shortdata(Ntraj,3,Xtrain));

% but use *all* data to get initial state (overwrite)
LDSparamsOBS.mu0 = mean(Xtrain)';
LDSparamsOBS.Info0 = inv(cov(Xtrain));

% test: filter, smoother
tic;
KFdstrbs    = KalmanFilter(LDSparamsOBS,gather(Ytest'));
SSerr       = sum((KFdstrbs.XHATMU' - Xtest).^2);
RsqKF       = 1 - SSerr./SStot;
RTSSdstrbs  = RTSsmoother(LDSparamsOBS,KFdstrbs);
SSerr       = sum((RTSSdstrbs.XHAT' - Xtest).^2);
RsqRTSS     = 1 - SSerr./SStot;
toc;

% plot?
if ~isempty(tvec)
    PosVelPlot(tvec,KFdstrbs.XHATMU',0);
    PosVelPlot(tvec,RTSSdstrbs.XHAT',0);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [LDSparams, ZhatF, ZhatS, RsqPF, RsqPS] = ...
    fullyObservedPoissonDecoder(trainData,testData,SStot)

keyboard

% train
LDSparams = learnfullyobservedLDS(trainData.Y,trainData.Z);
LDSparams = rmfield(LDSparams,'SigmaYX');
Nstates = size(LDSparams.A,1);
Nobsvs = size(LDSparams.C,1);
C = zeros(Nobsvs,1+Nstates,'like',trainData.Y);
for iObsv = 1:Nobsvs
    C(iObsv,:) = glmfitJGM(longdata(trainData.Z),...
        longdata(trainData.Y(:,iObsv,:)),'poisson');
end
muYX = C(:,1);
C = C(:,2:end);




% use data from *all* time to compute initial cumulants
xpctZ0 = mean(longdata(trainData.Z));
cvrnZ0 = cov(longdata(trainData.Z));

%%%% is this necessary?
LDSparams.muYX = muYX;
LDSparams.C = C;
LDSparams.mu0 = xpctZ0';
LDSparams.Info0 = inv(cvrnZ0);

% set the (particle) filter functions
Nparticles = 2500;
transitionFunc = @(ZZ,t)(LDSparams.A*ZZ + LDSparams.muX +...
    chol(LDSparams.SigmaX)'*randn(Nstates,Nparticles,'like',trainData.Y));
%%%% the repmat is unfortunate and can probably be eliminated---if
%%%% you're willing to rewrite poisspdf....
logEmissionFunc = @(YY,ZZ,t)(sum(poisslpdf(repmat(YY,[1,Nparticles]),...
    exp(C*ZZ + muYX)),1));


% test: filter, smoother
tic;
Ytest = squeeze(testData.Y(1,:,:));
Ztest = squeeze(testData.Z);
z0 = mvnrnd(xpctZ0,cvrnZ0,Nparticles)';
[Ztu,Wf] = particleFilter(Ytest,z0,transitionFunc,logEmissionFunc,1);
ZtuW = Ztu.*shiftdim(exp(Wf),-1);
ZhatF = permute(sum(ZtuW,2),[1,3,2]);
SSerr = sum((ZhatF - Ztest).^2,2)';
RsqPF = 1 - SSerr./SStot;
Ws = particleSmoother(Ztu,Wf,LDSparams,1);
ZtuW = Ztu.*shiftdim(exp(Ws),-1);
ZhatS = permute(sum(ZtuW,2),[1,3,2]);
SSerr = sum((ZhatS - Ztest).^2,2)';
RsqPS= 1 - SSerr./SStot;
toc;



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [LDSparams, ZhatF, ZhatS, RsqPF, RsqPS] = ...
    fullyObservedPolarDecoder(LDSparamsObs,trainData,testData,SStot)

keyboard

% replace the state-observation parameters
Ytrain = longdata(trainData.Y);
Ztrain = longdata(trainData.Z);
mvmtang = atan2(Ztrain(:,6),Ztrain(:,5));
speed = sqrt(sum(Ztrain(:,5:6).^2,2));
[C,RsqCV,Yres] = linregress([Ztrain(:,1:4),mvmtang,speed,Ztrain(:,7:9),...
    ones(size(speed,1),1)],Ytrain);
LDSparams = LDSparamsObs;
LDSparams.C = C(1:end-1,:)';
LDSparams.muYX = C(end,:)';
LDSparams.SigmaYX = cov(Yres);

% set the (particle) filter functions
Nparticles = 2500;
Nstates = size(LDSparams.A,1);
transitionFunc = @(ZZ,t)(LDSparams.A*ZZ + LDSparams.muX +...
    chol(LDSparams.SigmaX)'*randn(Nstates,Nparticles,'like',Ytrain));
logEmissionFunc = @(YY,ZZ,t)(mvnlpdf(repmat(YY,[1,size(ZZ,2)])',...
    (LDSparams.C*...
    [ZZ(1:4,:); atan2(ZZ(6,:),ZZ(5,:)); sqrt(sum(ZZ(5:6,:).^2).^2);ZZ(7:9,:)] +...
    LDSparams.muYX)',LDSparams.SigmaYX)');


% test: filter, smoother
tic;
Ytest = squeeze(testData.Y(1,:,:));
Ztest = squeeze(testData.Z(1,:,:));
Sigma0 = inv((LDSparams.Info0 + LDSparams.Info0')/2);
z0 = mvnrnd(LDSparams.mu0,Sigma0,Nparticles)';
[Ztu,Wf] = particleFilter(Ytest,z0,transitionFunc,logEmissionFunc,1);
ZtuW = Ztu.*shiftdim(exp(Wf),-1);
ZhatF = permute(sum(ZtuW,2),[1,3,2]);
SSerr = sum((ZhatF - Ztest).^2,2)';
RsqPF = 1 - SSerr./SStot;
Ws = particleSmoother(Ztu,Wf,LDSparams,1);
ZtuW = Ztu.*shiftdim(exp(Ws),-1);
ZhatS = permute(sum(ZtuW,2),[1,3,2]);
SSerr = sum((ZhatS - Ztest).^2,2)';
RsqPS = 1 - SSerr./SStot;
toc;



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [LDSparamsEM, Bxz_EM, LDSparamsObs, RsqStatic, RsqDynamic] =...
    latentStateLDSDecoder(Rtrain,Xtrain,Rtest,Xtest,params,...
    LDSparamsEM,SStot,tvec)


% learn a new system?
if isempty(LDSparamsEM)
    tic;
    Mstates = floor(params.numsUnits{1}/3);
    NepochsMax = 100;
    getTrajs = @(yrclass)(EFHdata2LDSdata(1,...
        @()(getSingleTrainingSequence('double',params))));
    LDSparamsEM = EM4LDS(Mstates,NepochsMax,'Normal',getTrajs,...
        'diagonal covariances',[1,1,0],...
        'parameter initialization','FactorAnalysis');
    toc
    
    % and now save it for future use
    if strcmp(params.datafile(1:6),'Chewie')
        filesuffix = [params.datafile(end-7:end-4),...
            params.datafile(end-11:end-8),'_01.mat'];
    else
        filesuffix = params.datafile(end-14:end);
    end
    
    subdir = 'scratch/';
    %%%subdir = '';
    save(sprintf('%sRBMish/BMI/%sLDSOrd%03d_%s_%s_%imsBins_%03dsTrainingtime_%s',...
        getdir('data'),subdir,Mstates,params.datatype,params.typeUnits{1}{1},...
        params.Nmsperbin,params.trainingtime,filesuffix),'LDSparamsEM');
end

% "observed variables"
LDSparamsEM.T   = size(Rtrain,1);
KFdstrbsTrainEM = KalmanFilter(LDSparamsEM,double(Rtrain'),'lightweight',1);
Vtrain          = [KFdstrbsTrainEM.XHATMU',ones(LDSparamsEM.T,1)];
clear KFdstrbsTrainEM Rtrain
LDSparamsEM.T   = size(Rtest,1);
KFdstrbsTestEM  = KalmanFilter(LDSparamsEM,double(Rtest'),'lightweight',1);
Vtest           = [KFdstrbsTestEM.XHATMU', ones(LDSparamsEM.T,1)];
clear KFdstrbsTestEM

% decode
[~,RsqStatic,XhatDynamic,RsqDynamic,Bxz_EM,LDSparamsObs] =...
    filterDecoder(Vtrain,Xtrain,Vtest,Xtest,1,SStot);

% plot?
if ~isempty(tvec), PosVelPlot(tvec,XhatDynamic,0); end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [RsqStatic,RsqDynamic] = linearREFHDecoder(...
    Rtrain,Xtrain,Qtrain,Rtest,Xtest,Qtest,wts,params,SStot,tvec)


% if necessary, train an rEFH
if isempty(wts)
    
    if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
    numsUnits = params.numsUnits;
    datagenargs = {};
    if isfield(params,'EACHBATCHISATRAJ')
        EACHBATCHISATRAJ = params.EACHBATCHISATRAJ;
    else
        EACHBATCHISATRAJ = 0;
    end
    if isfield(params,'Npretrain')
        Npretrain = params.Npretrain;
    else
        Npretrain = 0;
    end
    if isfield(params,'sparse')
        sparsitycost= params.sparse.cost;
        estimRate   = params.sparse.phidNewFrac;
        phidTarget  = params.sparse.phidTarget;
    else
        sparsitycost= 0;
    end
    getLatents  = params.getLatents;
    getData     = params.getData;
    NEFHs       = length(numsUnits)-1;
    Ncases      = params.Ncases;
    Nbatches    = params.Nbatches;
    Ntest       = params.NepochsMax + 1;
    NepochsMax  = params.NepochsMax;
    
    % init
    paramDisplay(params);
    wts = cell(NEFHs*2,1);
    
    % pretraining
    for iEFH = 1:NEFHs
        fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',...
            iEFH,sum(numsUnits{iEFH}),sum(numsUnits{iEFH+1}));
        RESTART = 1;
        
        % train
        tic; EFH; toc;
        
        % pack together weights for saving (hid => recog., vis => gener.)
        wts{iEFH} = [vishid; hidbiases];
        wts{NEFHs*2-iEFH+1} = [vishid'; visbiases'];
    end
    
    % and now save it for future use
    if strcmp(params.datafile(1:6),'Chewie')
        filesuffix = [params.datafile(end-7:end-4),...
            params.datafile(end-11:end-8),'_01.mat'];
    else
        filesuffix = params.datafile(end-14:end);
    end
    basedir = 'scratch';
    %%%basedir = 'spikecounts';
    save(sprintf('%sRBMish/EFHs/%s/wts_%s_%s_%iHid_%imsBins_%03dsTrainingtime_%s',...
        getdir('data'),basedir,params.datatype,params.typeUnits{1}{1},params.numsUnits{2},...
        params.Nmsperbin,params.trainingtime,filesuffix),'wts','params');
end

% test
[~,~,RsqStatic,XhatDynamic,RsqDynamic] = testEFHBMI(...
    Rtest,Xtest,Qtest,wts,params,Rtrain,Xtrain,Qtrain,SStot);

% plot?
if ~isempty(tvec), PosVelPlot(tvec,XhatDynamic,0); end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [UKFparams, RsqUKF] = UKFDecoder(Rtrain,Xtrain,T,Rtest,Xtest,...
    SStot,tvec,varargin)
% function [UKFparams, RsqUKF] = UKFDecoder(Rtrain,Xtrain,T,Rtest,Xtest,...
%   SStot,tvec,Ntapspast,Ntapspast);
%
% UKFdecoder wrapper file
%
% Rtrain - training spike rates  - Nsamples x Nneurons
% Xtrain - traning kinematics    - Nsamples x Nstates
% T      - num training samples  - Nsamples x 1
% Rtest  - testing spike rates   - Nsamples x Nneurons
% Xtest  - testing kinematics    - Nsamples x Nstates
% SStot  - sumsquares of testing - 1        x Nstates
%          SStot = sum((Xtest - mean(Xtest)).^2);
% tvec   - unused

%-------------------------------------------------------------------------%
% Revised: 04/18/17
%   -added arguments, code to control the tap order
% Created: 04/13/17
%   by JEO
%-------------------------------------------------------------------------%

% taps, etc.
fparams = [];  % extra params for f, not needed
ftaps = defaulter('ftaps',2,varargin{:});       % future taps
ptaps = defaulter('ptaps',3,varargin{:});       % past taps
mtaps = defaulter('mtaps',1,varargin{:});       % backbone; <=ftaps+ptas
htaps = defaulter('htaps',1,varargin{:});       % spiking history taps
alphamax = 16; % maximum regularization (ridge regression) parameter

% f: backbone (px py vx vy ax ay) -> variables encoded linearly by spikes
f = @(a,b)([...
    a;...
    sqrt(a(1,:).^2+a(2,:).^2);...
    sqrt(a(3,:).^2+a(4,:).^2);...
    sqrt(a(5,:).^2+a(6,:).^2);...
    a(1,:).*a(3,:);...
    a(2,:).*a(4,:);...
    ]);

% fit
params = fit_ar_ukf_hist_rrcv(Xtrain,Rtrain,f,fparams,...
    ftaps,ptaps,mtaps,htaps,alphamax);

% re-fit initial states from *all* training data
xinit = mean(Xtrain);
vinit = cov(Xtrain);

% test
[XhatUKF, ~] = ar_ukf_hist(Rtest,params,xinit,vinit,1);

% compute R^2 using predictions for *current* kinematic vars
Nstates = size(Xtrain,2);
indsNow = (1:Nstates) + (ptaps-1)*Nstates;
RsqUKF = 1 - sum((XhatUKF(:,indsNow)- Xtest).^2)./SStot;

% todo fill these
UKFparams.A = params.F;
UKFparams.SigmaX = [];
UKFparams.nuX = [];
UKFparams.C = [];
UKFparams.SigmaYX = [];
UKFparams.muYX = [];
UKFparams.mu0 = [];
UKFparams.Cvrn0 = [];
UKFparams.T = [];
UKFparams.params = params;

if ~isempty(tvec), PosVelPlot(tvec,XhatUKF,0); end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function RsqsAvg = averageDecentRsqs(SSerr,SStot,fignum)

Rsqs = 1 - SSerr./SStot;
%figure(fignum)
%imagesc(reshape(Rsqs,[4,292/4]));
%colorbar

foo = sort(Rsqs);
RsqsAvg = mean(foo(50:end));


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Rtrain,Xtrain,Qtrain] = getSingleTrainingSequence(dataclass,params)

[Xtrain,Qtrain] = params.getLatents([],dataclass,'sequencelength','singlesequence');
[Rtrain,Qtrain] = params.getData(Xtrain,Qtrain);

end
%-------------------------------------------------------------------------%

% NOTES/TO DO
% (1) Merge all the files that do KFs on neural data, with switch
%   statements for the different preprocessing
%   --at least clean up all the unused garbage in KF4HHS
%   --it also seems like you should be able to eliminate some of the HHS
%   fxns in favor of some of the ../tools for Kalman filters....
% (2) Consider moving from pos/vel/acc to pos(t)/pos(t+1)/pos(t+2)

