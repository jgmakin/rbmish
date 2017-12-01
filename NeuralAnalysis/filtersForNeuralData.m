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
Nsplot = 60;
Nmsperbin = params.Nmsperbin;
Nmspers = 1000; % fact
if TOPLOT, tvec = 1:ceil(Nsplot/(Nmsperbin/Nmspers)); else, tvec=[]; end

% get data
[Rtrain,Xtrain,Qtrain] = getSingleTrainingSequence(dataclass,params);
[Rtest,Xtest,Qtest] = params.getTestData(dataclass);
SStot = sum((Xtest - mean(Xtest)).^2);
NdataTest = size(Rtest,1);
if ~isempty(tvec), PosVelPlot(tvec,Xtest,1); end


% (1) (static) linear regression
if any(strcmp(decodernames,'static'))
    [~,~,RsqStatic] = staticLinearDecoder(Rtrain,Xtrain,...
        Rtest,Xtest,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'static'))       = RsqStatic;
end


% (2) observed LDS filtered with KF
if any(strcmp(decodernames,'kfobs'))
    [LDSparamsObs,~,RsqKF_Obs] = ...
        fullyObservedLDSDecoder(Rtrain,Xtrain,Qtrain.T,Rtest,Xtest,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'kfobs'))        =...
        RsqKF_Obs(:,:,:,any(strcmp(decodernames,'kfobs')));
    %Rsqs(:,:,:,strcmp(decodernames,'rtssobs'))      =...
    %   RsqRTSS_Obs(:,:,:,any(strcmp(decodernames,'rtssobs')));
end


% (3,4) KF with latent state + (static map | KF)
if any(ismember(decodernames,{'kfemstatic','kfemdynamic'}))
    
    NhiddenStates = floor(params.numsUnits{1}/3);
    
    % training-data generator function
    getTrainingData = @(yrclass)(EFHdata2LDSdata(1,...
        @()( deal(Rtrain,Xtrain) )));
    [LDSparamsEM,Bxz_EM,LDSparamsObs,~,RsqKF_EMstatic,~,RsqKF_EMdynamic] =....
        latentStateLDSDecoder(getTrainingData,Rtrain,Xtrain,Rtest,Xtest,...
        params.Nmsperbin,NhiddenStates,LDSparamsEM,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'kfemstatic'))   =...
        RsqKF_EMstatic(:,:,:,any(strcmp(decodernames,'kfemstatic')));
    Rsqs(:,:,:,strcmp(decodernames,'kfemdynamic'))  =...
        RsqKF_EMdynamic(:,:,:,any(strcmp(decodernames,'kfemdynamic')));

    % and now save it for future use
    if strcmp(params.datafile(1:6),'Chewie')
        filesuffix = [params.datafile(end-7:end-4),...
            params.datafile(end-11:end-8),'_01.mat'];
    else
        filesuffix = params.datafile(end-14:end);
    end
    
    subdir = 'scratch/';
    %%%subdir = 'EMparams/';
    save(sprintf('%sRBMish/BMI/%sLDSOrd%03d_%s_%s_%imsBins_%03dsTrainingtime_%s',...
        getdir('data'),subdir,NhiddenStates,params.datatype,...
        params.typeUnits{1}{1},params.Nmsperbin,params.trainingtime,...
        filesuffix),'LDSparamsEM','LDSparamsObs','Bxz_EM');
end


% (5) UKF
if any(strcmp(decodernames,'ukf'))
    [UKFparams,~,RsqUKF] = UKFDecoder(Rtrain,Xtrain,Qtrain.T,...
        Rtest,Xtest,SStot,'ftaps',0,'ptaps',1,'htaps',0);
    Rsqs(:,:,:,strcmp(decodernames,'ukf'))          = RsqUKF;
end


% (6) UKF with taps
if any(strcmp(decodernames,'ukfwithtaps'))
    ptaps = floor(128/Nmsperbin) + 1;
    ftaps = floor(128/Nmsperbin);
    [UKFparams, RsqUKFwithtaps] = UKFDecoder(Rtrain,Xtrain,Qtrain.T,...
        Rtest,Xtest,SStot,'ftaps',ftaps,'ptaps',ptaps);
    Rsqs(:,:,:,strcmp(decodernames,'ukfwithtaps'))  = RsqUKFwithtaps;
end


% (7) Wiener filter
if any(strcmp(decodernames,'wf'))
    [betaWF, RsqWF] = WFDecoder(Rtrain,Xtrain,Qtrain.T,params.Nmsperbin,...
        Rtest,Xtest);
    % NB that we don't pass SStot to this decoder because the testing set
    % will be a few samples (Nbin_of_history) shy of the full set.
    Rsqs(:,:,:,strcmp(decodernames,'wf'))          = RsqWF;
end


% (8,9) rEFH
if any(ismember(decodernames,{'refhstatic','refhdynamic',...
        'refhstatic_stdnrml','refhdynamic_stdnrml'}))
    
    [wts,~,RsqREFHstatic,~,RsqREFHdynamic] = linearREFHDecoder(Rtrain,Xtrain,...
        Qtrain,Rtest,Xtest,Qtest,wts,params,SStot);
    
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
    
    Rsqs(:,:,:,strcmp(decodernames,'refhstatic'))   =... 
        RsqREFHstatic(:,:,:,any(strcmp(decodernames,'refhstatic')));
    Rsqs(:,:,:,strcmp(decodernames,'refhdynamic'))  =... 
        RsqREFHdynamic(:,:,:,any(strcmp(decodernames,'refhdynamic')));
    Rsqs(:,:,:,strcmp(decodernames,'refhstatic_stdnrml'))   =... 
        RsqREFHstatic(:,:,:,any(strcmp(decodernames,'refhstatic_stdnrml')));
    Rsqs(:,:,:,strcmp(decodernames,'refhdynamic_stdnrml'))  =... 
        RsqREFHdynamic(:,:,:,any(strcmp(decodernames,'refhdynamic_stdnrml')));
end



% (10,11) Poisson emissions (particle filter)
if any(ismember(decodernames,{'particlefilter','particlesmoother'}))
    [LDSparamsPE, PFdstrbs, PSdstrbs, RsqPF, RsqPS] = ...
        fullyObservedPoissonDecoder(trainData,testData,SStot);
    Rsqs(:,:,:,strcmp(decodernames,'particlefilter'))   =... 
        RsqPF(:,:,:,any(strcmp(decodernames,'particlefilter')));
    Rsqs(:,:,:,strcmp(decodernames,'particlesmoother')) =... 
        RsqPS(:,:,:,any(strcmp(decodernames,'particlesmoother')));
end

% (12,13) observations linear in polar transform of vel (particle filter)
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
    Ytest,Xtest,SStot)

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

