function eStats = mastertest(wts,params,varargin)
% Compute, in the simplest way possible, the optimal posterior mean for the
% various models.
%
% USAGES:
%{

    % '2Dinteg'   
    load([getdir('data'),'RBMish/EFHs/new/wts_2Dinteg_900.mat'])
    eStats = mastertest(wts,params);

    % '2Dinteg' but with gains covering the whole range
    load([getdir('data'),'RBMish/EFHs/new/wts_2Dinteg_AllGains_900.mat'])
    params.swing = 0.8; % don't *test* with gains = 0.
    eStats = mastertest(wts,params);

    % '1Dinteg'
    load([getdir('data'),'RBMish/EFHs/new/wts_1Dinteg_60.mat'])
    eStats = mastertest(wts,params);

    % '1Daddition'  
    load([getdir('data'),'RBMish/EFHs/new/wts_1Daddition_160.mat'])
    eStats = mastertest(wts,params,'posteriorlist',...
        {'unisensory','optimal','EFHuni','EFHmulti'})

    % 'HierL2' (hierarchical integration)
    load([getdir('data'),'RBMish/EFHs/new/wts_HierL2_900.mat'])
    params = setParams('datatype','HierL2');
    eStats = mastertest(wts,params);

    % '2Dtwoarms'
    load([getdir('data'),'RBMish/EFHs/new/wts_2Dtwoarms_900.mat'])
    eStats = mastertest(wts,params);

    % 'rEFH_1D_LTI-PPC'
    load([getdir('data'),'RBMish/EFHs/new/wts_rEFH_1D_LTI-PPC_240.mat'])
    yrposteriors = {'unisensory','EM1','EFHuni','EM2','optimal'};
    eStats = mastertest(wts,params,'LDS_EM',1,'LDS_EM',2,...
        'posteriorlist',yrposteriors)
        
    % 'TRBM_1D_LTI-PPC'
    load([getdir('data'),'\RBMish\EFHs\new\wts_TRBM_1D_LTI-PPC_240.mat'])
    yrposteriors = {'unisensory','EM1','EFHuni','EM2','optimal'};
    eStats = mastertest(bestwts,params,'LDS_EM',1,'LDS_EM',2,...
        'posteriorlist',yrposteriors)

    % 'RTRBM_1D_LTI-PPC'
    load([getdir('data'),'\RBMish\EFHs\new\wts_RTRBM_1D_LTI-PPC_240.mat'])
    yrposteriors = {'unisensory','EM1','EFHuni','EM2','optimal'};
    eStats = mastertest(wts,params,'LDS_EM',1,'LDS_EM',2,...
        'posteriorlist',yrposteriors)

    % '1DrEFHwithEC'
    load([getdir('data'),'RBMish/EFHs/new/wts_rEFH_LTI-PPC_with_EC_480.mat'])
    yrposteriors = {'unisensory','EM2','EFHuni','EM3','optimal'};
    eStats = mastertest(wts,params,'LDS_EM',2,'LDS_EM',3,...
        'posteriorlist',yrposteriors)





    % '2Dinteg' but with a non-flat prior distribution of stimuli
    load([getdir('data'),'RBMish/EFHs/wts_2Dinteg_NonFlatPrior_900_28-Dec-2016.mat'])
    eStats = mastertest(wts,params)

    % '2DintegDecoupled' (occasionally decoupled 2D integration)
 	load([getdir('data'),'RBMish/EFHs/wtsMixture120524.mat'])
    eStats = mastertest(wts,params,'correlation',0.7,'posteriorlist',...
        {'unisensory','optimal','EFHmulti'})

   
%}

%-------------------------------------------------------------------------%
% Revised: 07/09/14
%   -made it work for params.MODEL = '2DintegDecoupled'
% Revised: 07/07/14
%   -made it work for params.MODEL = 'HierL2'
% Revised: 07/05/14
%   -made it work for params.MODEL = '1Daddition'
% Revised: 07/04/14
%   -made it work, particulary with tensor operations (no loops)
%   -made it work for params.MODEL = '1Dinteg'
% Created: 06/30/14
%   by JGM
%-------------------------------------------------------------------------%

%%%% TO DO
% (3) allow the user to request:
%       (b) which line styles, also (cf. HierL2, and maybe coord trans)?
%       (c) the filename for the figure to be saved under?
% (4) Create a mastertrain.m?  That would either replace or call DBNtrain
%   and recurrentEFHtrain.
% (5) allow *multiple* weights?  So e.g. to compare RTRBM, TRBM, and rEFH?
%%%%%



% params
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[updownargs, whichposteriors, extraDecoders, params] = setAudibles(params,varargin{:});


% make data
tic
[R0,X,Q] = params.getTestData(dataclass);
toc

% calculate posterior distributions for original and EFHed data
tic
fprintf('\ncalculating posterior distributions...\n');
[unisensCmlnts0,unisensCmlnts1,multisensCmlnts] =...
    getPosteriorCumulants(R0,X,Q,wts,params,updownargs);
altCmlnts = altDecoders(unisensCmlnts0,unisensCmlnts1,Q,extraDecoders,params);
toc

% collect the useful posteriors, compute errors, and display
tic
fprintf('\ncalculating error statistics...\n')
testPosteriors = assembleTestPosteriors(unisensCmlnts0,unisensCmlnts1,...
    multisensCmlnts,altCmlnts,params,whichposteriors{:});

Sneutral = Q.latent2stim{strcmp(params.mods,params.NS)}(X);
eStats = getErrorStats(testPosteriors,Sneutral);
dispErrStats(eStats,params.NS);
fprintf('\ndone\n');
toc




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [updownargs, whichposteriors, decoders, params] =...
    setAudibles(params,varargin)

% defaults
if isfield(params,'dynamics')
    updownargs = defaulter('propagation',{'means'},varargin{:});
    %%% that's how it is in the LtTMSwPC
else
    updownargs = defaulter('propagation',{'Nsamples'},varargin{:});
    %%% that's how it is in the MIvDE
end
params.smpls = defaulter('numsamples',15,varargin{:});
whichposteriors = defaulter('posteriorlist',...
    {'unisensory','optimal','EFHuni'},varargin{:});

% "audibles"
datagenargs = {};   iD = 1;
decoders    = {};   iE = 1;
for iArg = 1:2:length(varargin)
    switch varargin{iArg}
        case {'neuralNetwork','LDS_EM'}
            decoders{iE} = varargin{iArg};
            decoders{iE+1} = varargin{iArg+1};
            iE=iE+2;
        otherwise
            %%% do nothing
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [unisensCmlnts0,unisensCmlnts1,multisensCmlnts] =...
    getPosteriorCumulants(R0,X,Q,wts,params,updownargs)

% decode; compute optimal (multisensory); store
[Shat0,ttlSpks0] = decodeDataPPC(R0,X,Q,params);
Info0 = GTPNposteriorInfo(ttlSpks0,params);
unisensCmlnts0 = cumulantNeutralize(Shat0,Info0,params);

% benchmarks: unisensory and multisensory (integrated) estimates
if isfield(params,'dynamics')
    LDSparamsTrue = getLDSparams(params.dynamics);
    multisensCmlnts = KFposteriorization(unisensCmlnts0,Q,LDSparamsTrue,params);
    
    % store
    %%% Consider making this unnecessary.  But there's no great option.
    %%% You'd like to return posteriors over both Joint-Angle and
    %%% Efference-Copy (where applicable), but mastertest only produces
    %%% plots for params.NS. One radical plan would be to make params.NS
    %%% hold *both* of these.... For now, a hack.  Some of this might also
    %%% be fixable by putting in NaNs, and then relying on the plotting
    %%% fxns to ignore nan stats.....  But note that things work right for
    %%% LDStest.m.
    NSind = strcmp(unisensCmlnts0.srcs,params.NS);
    multisensCmlnts.Xpct = multisensCmlnts.Xpct(:,:,NSind);
    multisensCmlnts.Info = multisensCmlnts.Info(:,:,:,NSind);
else
    if strcmp(params.datatype,'HierL2')
        [Shat0,ttlSpks0] = decodeDataPPC(Q.R1and2,X,Q.Q0,Q.params0);
        Info0 = GTPNposteriorInfo(ttlSpks0,Q.params0);
        unisensCmlntsL1 = cumulantNeutralize(Shat0,Info0,Q.params0);
        unisensCmlnts0.Xpct = cat(3,unisensCmlnts0.Xpct,unisensCmlntsL1.Xpct);
        unisensCmlnts0.Info = cat(4,unisensCmlnts0.Info,unisensCmlntsL1.Info);
        unisensCmlnts0.srcs = [unisensCmlnts0.srcs,unisensCmlntsL1.srcs];       
    end
    
    % decode; compute optimal (multisensory); store
    multisensCmlnts = gaussPosteriorization(unisensCmlnts0);
end

% decode the "updated" input layer
[~,~,Shat1,Info1] = testEFHPPC(R0,X,Q,wts,params,'updownargs',updownargs{:});
unisensCmlnts1 = cumulantNeutralize(Shat1,Info1,params);




end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function altCmlnts = altDecoders(unisensCmlnts0,unisensCmlnts1,Q,decoders,params)

% initialize Xpct
if isempty(decoders)
    altCmlnts.Xpct = [];
else
    altCmlnts.Xpct = zeros(size(unisensCmlnts1.Xpct,1),params.Ndims,length(decoders)/2);
    altCmlnts.src = cell(1,length(decoders)/2);
    for iD = 1:2:length(decoders)
        switch decoders{iD}
            
            %%%%% BROKEN
            case 'neuralNetwork'
            %%%%
                
            case 'LDS_EM'
                                
                [altCmlnts.Xpct(:,:,(iD+1)/2),altCmlnts.src{(iD+1)/end}] =...
                    filterWithSavedEMbasedModel(unisensCmlnts1,Q,...
                    decoders{iD+1},...
                    any(strcmp(params.mods,'Efference-Copy')),...
                    params.Ndims,params.datatype);
                
            case 'EFHmulti'
                [altCmlnts.Xpct(:,:,(iD+1)/2),altCmlnts.src{(iD+1)/end}] =...
                    gaussPosteriorization(unisensCmlnts1);
                
            otherwise
                error('unrecognized decoder -- jgm\n');
        end
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function testPosteriors = assembleTestPosteriors(unisensCmlnts0,...
    unisensCmlnts1,multisensCmlnts,altCmlnts,params,varargin)
% Assemble the posteriors to be tested and plotted!

% "audibles"
iPost = 0;
jPost = 0;
for iArg = 1:length(varargin)
    switch varargin{iArg}
        case 'unisensory'
            for iMod = 1:size(unisensCmlnts0.Xpct,3)
                %%% HACK!
                if strcmp(unisensCmlnts0.srcs{iMod},'Efference-Copy')
                else
                iPost = iPost + 1;
                testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,iMod);
                testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{iMod};
                end
            end
            
        case 'Joint-Angle'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,...
                strcmp(params.mods,'Joint-Angle'));
            testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{...
                strcmp(params.mods,'Joint-Angle')};
        case 'Hand-Position'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,...
                strcmp(params.mods,'Hand-Position'));
            testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{...
                strcmp(params.mods,'Hand-Position')};
        case 'optimal'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts.Xpct;
            testPosteriors.srcs{iPost} = 'opt';
        case 'EFHuni'
            iPost = iPost + 1;
            indsN = strcmp(params.NS,unisensCmlnts1.srcs);
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts1.Xpct(:,:,indsN);
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHmulti'
            iPost = iPost + 1;
            fprintf('\ndecoding EFH by integrating "updated" populations\n');
            multisensCmlnts1 = gaussPosteriorization(unisensCmlnts1);
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts1.Xpct;
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHnn'
            iPost = iPost + 1;
            jPost = jPost + 1;
            fprintf('\ndecoding EFH with ANN\n');
            testPosteriors.Xpct(:,:,iPost) = altCmlnts.Xpct(:,:,jPost);
            testPosteriors.srcs{iPost} = 'EFHnn';
        case {'EM1','EM2','EM3'}
            iPost = iPost + 1;
            jPost = jPost + 1;
            testPosteriors.Xpct(:,:,iPost) = altCmlnts.Xpct(:,:,jPost);
            testPosteriors.srcs{iPost} = altCmlnts.src{jPost};
        otherwise
            error('undefined posteriors! -- jgm\n');
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Xpct,name] = filterWithSavedEMbasedModel(unisensCmlnts0,Q,...
    NKFstates,WITH_EC,Ndims,datatype)
% load EM-derived model and get its estimates

% assemble file name
fileprefix = [getdir('data'),'RBMish/'];
if WITH_EC
    modelsuffix = '_withEC';
else
    modelsuffix = '';
end

filename = [fileprefix,'EMparams/LDSOrd',num2str(NKFstates),'_',...
    num2str(Ndims),'D','_',datatype,modelsuffix,'_VaryingHids.mat'];
load(filename,'allparams','params');

% EM or EM best?
Nstates = size(params.dynamics.C,2);
if Nstates>NKFstates, name='EM'; else name='EM (best)'; end

% Kalman filter
multisensCmlnts = KFposteriorization(unisensCmlnts0,Q,allparams(1),params);
NSind = strcmp(unisensCmlnts0.srcs,params.NS);
Xpct = multisensCmlnts.Xpct(:,:,NSind);

end
%-------------------------------------------------------------------------%
