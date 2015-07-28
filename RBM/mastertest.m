function eStats = mastertest(wts,params,varargin)
% Compute, in the simplest way possible, the optimal posterior mean for the
% various models
%
% USAGES:
%{

  % '2Dinteg'
  load('results\finalwts\wtsStandard140613.mat','wts','params');
  eStats = mastertest(wts,params);

  % '1Dinteg'
  load('results\finalwts\wts1Dinteg140704.mat','wts','params');
  eStats = mastertest(wts,params);

  % '1Daddition'
  load('results\finalwts\wtsCoordTrans140616.mat','wts','params');
  eStats = mastertest(wts,params,'posteriorlist',...
      {'unisensory','optimal','EFHmulti'})

  % 2D integration with prior distribution
  load('results\finalwts\wts2dIntegWithPrior131216.mat','wts','params','p0');
  params.p0 = p0;
  mastertest(wts,params)

  % 'HierL2' (hierarchical integration)
  load('results\finalwts\wtsStandard140613','wts','params');
  params.smpls = 15; wts0 = wts; params0 = params;
  load('results\finalwts\wtsStandardL2140616','wts','params');
  mastertest(wts,params,'prevweights',wts,'prevparams',params)

  % '2DintegDecoupled' (occasionally decoupled 2D integration)
  load('results\finalwts\wtsMixture120524.mat')
  mastertest(wts,params,'correlation',0.7,'posteriorlist',...
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
% (2) allow the user to request:
%       (b) which line styles, also (cf. HierL2, and maybe coord trans)?
%       (c) the filename for the figure to be saved under?
% (3) Create a mastertrain.m 
% (4) consider fixing all the other m-files that generated figures/useful
% results for your paper, making them use the new formats and code.  Then
% you may be able to retire a lot of files (see fakah.m)
%%%%%

% params
[~,machine] = system('hostname');
params.machine = strtrim(machine);
[updownargs, datagenargs, whichposteriors, NNdecoder, params] =...
    setAudibles(params,varargin{:});

% make data
tic
Nbatches = 40000/params.Ncases;
[R0,S,~,Q] = DATAGENPP(Nbatches,params,datagenargs{:});
[R0,S] = longdata(R0,S);
toc

%%%%%%
% load('results\DARPA\deadinds.mat','deadIndsOneQrtr','deadIndsOneHalf',...
%     'deadIndsThreeQrtrs','deadIndsFourQrtrs'); 
% deadInds = {deadIndsOneQrtr,deadIndsOneHalf,deadIndsThreeQrtrs,...
%     deadIndsFourQrtrs};
% R0(:,deadInds{4}) = 0;
%%%%%%


% send up and down through the EFH
tic
fprintf('\nintegrating with the EFH...\n');
[~,R1] = updown(R0,wts,params,updownargs{:});
fprintf('\nmaking any special-case adjustments to the data...\n');
[R0,R1,S,params] = dataAdjust(R0,R1,S,Q,params);
toc
%%% speed this up (is there a loop in here?)


% calculate posterior distributions for original and EFHed data
tic
fprintf('\ncalculating posterior distributions...\n');
[unisensCmlnts0,multisensCmlnts0] = getPosteriorCumulants(R0,params);
[unisensCmlnts1,multisensCmlnts1] = getPosteriorCumulants(R1,params);
nnXpct = getNNxpct(R0,R1,NNdecoder,params);
toc

% collect the useful posteriors, compute errors, and display
tic
fprintf('\ncalculating error statistics...\n')
testPosteriors = assembleTestPosteriors(unisensCmlnts0,unisensCmlnts1,...
    multisensCmlnts0,multisensCmlnts1,nnXpct,params,whichposteriors{:});
eStats = getErrorStats(testPosteriors,S(:,:,strcmp(params.mods,params.NS)));
dispErrStats(eStats,params.NS);
fprintf('\ndone\n');
toc




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [updownargs, datagenargs, whichposteriors, NNdecoder, params] =...
    setAudibles(params,varargin)

% defaults
updownargs = {'Nsamples'};
params.smpls = 15;
datagenargs = {}; iD = 1;
whichposteriors = {'unisensory','optimal','EFHuni'};
NNdecoder = {};

% "audibles"
for iArg = 1:2:length(varargin)
    switch varargin{iArg}
        case 'propagation'
            updownargs{1} = varargin{iArg+1};
        case 'numsamples'
            params.smpls = varargin{iArg+1};
        case 'posteriorlist'
            whichposteriors = varargin{iArg+1};
        case 'neuralNetwork'
            NNdecoder = varargin{iArg+1}; 
        otherwise
            datagenargs{iD} = varargin{iArg};
            datagenargs{iD+1} = varargin{iArg+1};
            iD=iD+2;
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [R0,R1,S,params] = dataAdjust(R0,R1,S,Q,params)

switch params.MODEL
    
    case 'HierL2'
        % put all the input populations together
        params.mods = [Q.paramsOLD.mods params.mods{:}];
        
        Nvis = params.numsUnits(1);
        switch params.typeUnits{1}
            case 'PB', inds3 = 1:(Nvis-params.t);
            case 'BP', inds3 = (params.t+1):Nvis;
            otherwise, error('unexpected unit types in hierarchy -- jgm');
        end
        R0 = cat(2,longdata(Q.D),R0(:,inds3,:));
        R1 = cat(2,longdata(Q.D),R1(:,inds3,:));
        
        S = cat(3,longdata(Q.S),S);
        
        params.smin = cat(2,Q.paramsOLD.smin,params.smin);
        params.smax = cat(2,Q.paramsOLD.smax,params.smax);
        params.Nmods = Q.paramsOLD.Nmods + params.Nmods;
        params.posmin = Q.paramsOLD.posmin;
        params.posmax = Q.paramsOLD.posmax;
        %%% change params.numsUnits??
        
    otherwise
        fprintf('\nno adjustments necessary...\n')
end
end
%-------------------------------------------------------------------------%  

%-------------------------------------------------------------------------%
function [unisensoryCumulants,multisensoryCumulants] =...
    getPosteriorCumulants(R,params)

% kind of ugly
if strcmp(params.MODEL,'2DintegDecoupled')
    DECOUPLED = (R(:,params.N) > 0.08);
    R(:,params.N) = 0;
else
    DECOUPLED = false(size(R,1),1);
end

[cntrOfMass, ttlSpks] = GTPNsuffstats(R,params);
ttlSpks(DECOUPLED,~strcmp(params.NS,params.mods)) = 0;              %%
Info = GTPNposteriorInfo(ttlSpks,params);
%%%%%%
fprintf('\n\ndoing terrible clamping-at-edges thing to avoid complex numbers\n\n');
smin = shiftdim(params.smin,-1);
smax = shiftdim(params.smax,-1);
cntrOfMass = bsxfun(@ge,cntrOfMass,smin).*cntrOfMass +...
    bsxfun(@times,bsxfun(@lt,cntrOfMass,smin),smin);
cntrOfMass = bsxfun(@le,cntrOfMass,smax).*cntrOfMass +...
    bsxfun(@times,bsxfun(@gt,cntrOfMass,smax),smax);
%%%%%%
unisensoryCumulants = cumulantNeutralize(cntrOfMass,Info,params);
unisensoryCumulants.Xpct(DECOUPLED,:,...
    ~strcmp(unisensoryCumulants.srcs,params.NS)) = NaN;             %%
multisensoryCumulants = gaussPosteriorization(unisensoryCumulants);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function nnXpct = getNNxpct(R0,R1,NNdecoder,params)

if ~isempty(NNdecoder)
    switch NNdecoder.userdata
        case 'R1toNS'
            %%% bad reduplication...
            [cntrOfMass, ttlSpks] = GTPNsuffstats(R1,params);
            [Nexamples,Ndims,Nmods] = size(cntrOfMass);
            X = [reshape(cntrOfMass,[Nexamples,Ndims*Nmods]),ttlSpks]';
            nnXpct = NNdecoder(X)';
        case 'V0toVel'
            % ...
        case 'R0toSomething'
            % ...
        otherwise
            error('unrecognized neural-network decoder -- jgm\n');
    end    
else
    nnXpct = [];
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function testPosteriors = assembleTestPosteriors(unisensCmlnts0,...
    unisensCmlnts1,multisensCmlnts0,multisensCmlnts1,nnXpct,params,varargin)
% Assemble the posteriors to be tested and plotted!


% "audibles"
iPost = 0;
for iArg = 1:length(varargin)
    switch varargin{iArg}
        case 'unisensory'
            for iMod = 1:size(unisensCmlnts0.Xpct,3)
                iPost = iPost + 1;
                testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,iMod);
                testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{iMod};
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
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts0.Xpct;
            testPosteriors.srcs{iPost} = 'opt';
        case 'EFHuni'
            iPost = iPost + 1;
            indsN = strcmp(params.NS,unisensCmlnts1.srcs);
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts1.Xpct(:,:,indsN);
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHmulti'
            iPost = iPost + 1;
            fprintf('\ndecoding EFH by integrating "updated" populations\n');
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts1.Xpct;
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHnn'
            iPost = iPost + 1;
            fprintf('\ndecoding EFH with ANN\n');
            testPosteriors.Xpct(:,:,iPost) = nnXpct;
            testPosteriors.srcs{iPost} = 'EFH';
        otherwise
            error('undefined posteriors! -- jgm\n');
    end
end


end
%-------------------------------------------------------------------------%
