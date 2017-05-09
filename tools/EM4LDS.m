function qLDSparams = EM4LDS(Nstates,NepochsMax,datadstrb,getTrajs,varargin)
% Expectation-maximization for an LTI system
%
% USAGES:
%
%   LDSparams = EM4LDS(Nstates,NepochsMax,'Normal',...
%       @(yrclass)(EFHdata2LDSdata(Ntraj,T,yrclass,params)));
%
%   LDSparams = EM4LDS(Nstates,NepochsMax,'Poisson',...
%       @(yrclass)(EFHdata2LDSdata(Ntraj,T,yrclass,params)),...
%       'diagonal covariances',[1,1,0],...
%       'inferencetype','particle',...
%       'state initialization',getInitialStateWrapper);
%
% EM4LDS(Nstates,params) fits a linear, time-invariant, dynamical system
% with expectation maximization (EM).  Training data are generated with the
% function passed in as the fourth argument.  Inference is exact (assuming
% jointly Gaussian data) unless particle methods are requested.  In this 
% case, the emission can (also) be 'GammaFixedScale' or 'Poisson'.


%-------------------------------------------------------------------------%
% Revised: 09/09/16
%   -added case 'GammaFixedScale' for particle filtering.
%   -changed what data are collected in the particle E step (to be handed
%   to the new particleIRLS.m)
% Revised: 09/02/16
%   -incorporated particle filtering and smoothing as an alternative
%   version of the E step (formerly in particleEM4LDS.m).
%   -made other minor changes to accommodate this merge
% Revised: 08/17/15
%   -large revision: rewrote to accommodate LDSdata that has
%   variable-length trajectories (and hence is an array structure rather
%   than a structure of tensors)
% Revised: 08/26/14
%   -eliminated "noisily observed controls"
% Revised: 04/01/14
%   -forced the inference algorithms to use prestored SigmaV & SigmaYX---if
%   they exist
%   -added plotting function for log-likelihoods
%   -added log likelihood (negative cross entropy) for p(v)
% Revised: 03/10/14
%   -made some adjustments for 1D data
% Revised: 01/10/14
%   -rewrote initLDSparams to accommodate controls; fixed some of its bugs
% Revised: 01/07/14
%   -eliminated outputs eKF and eRTSS
% Revised: 12/02/13
%   -changed EM loop to hold off update until after all Ncases
%   trajectories; then added outer loop through "epochs"
% Revised: 11/25/13
%   -adapted gatherXpctSuffStats to allow for nonzero residual means in the
%   emission (something like the column of ones...)
% Revised: 10/30/13
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%




%%%% TO DO
% (0) For the particle filter, new data (at Nrenew) can drastically
% increase the cross entropy.  (Your temporory solution is not to renew the
% data.)  How now?  
% (1) Here's a problem: for the emission types that require IRLS, and
% therefore require learning on ONE TRAJECTORY AT A TIME, estimates of the
% initial position converge quickly to have very small covariances---which
% is wrong.  If fresh data are then created (e.g. at NepochsRenew epochs), 
% the inference algorithm will insist (wrongly) that the new trajectory
% starts at the old one's starting location.
%   -One possibility is to over-write x0 and cov0 for these cases, using
%   instead *all* the data.



% numbers
NepochsRenew = 1; %%% consider abolishing
thr = 0.0005;	% 1/20 of a percent

% variable arguments
iParams     = find(strcmp(varargin,'params')); %%% confusing; should rename
iData       = find(strcmp(varargin,'data'));
DIAGCOVS    = defaulter('diagonal covariances',[0,0,0],varargin{:});
USELOGS     = defaulter('uselogs',1,varargin{:});
INFRTYPE    = defaulter('inferencetype','exact',varargin{:});
TOPLOT      = defaulter('toplot',0,varargin{:});
pLDSparams  = defaulter('true dynamics',[],varargin{:});
initscheme  = defaulter('parameter initialization','random',varargin{:});
stateinitfunc = defaulter('state initialization',{},varargin{:});

% for particle-based approaches
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
if strcmp(INFRTYPE,'particle')
    if ~isempty(stateinitfunc)
        error('Need a state initialization function for particle filter');
    end
    
    Nparticles = 2000; % 500;
    if Nparticles > 40000
        HADBEENCLOSED = isempty(gcp('nocreate'));
        if HADBEENCLOSED, parpool; end
    end
else
    Nparticles = 0;
end

% EXPECTATION-MAXIMIZATION
for iEpoch = 1:NepochsMax
    
    %------------------- DATA GENERATION/INITIALIZATION ------------------%
    % every Nrenew epochs...
    if mod(iEpoch-1,NepochsRenew)==0
        
        if isempty(iData)
            trajs = getTrajs(dataclass);
        else
            trajs = varargin{iData+1};
        end
        
        % initialize model params
        if iEpoch==1
            if isempty(iParams)
                qLDSparams = initLDSparams(trajs,Nstates,datadstrb,...
                    initscheme,DIAGCOVS);
            else
                qLDSparams = varargin{iParams+1};
            end
            XNtrpYold = Inf;
            
            % malloc
            XNtrpY = zeros(1,NepochsMax,'like',trajs(1).Y);
            
            % figure preparation
            TOPLOT = TOPLOT&usejava('desktop');
            if TOPLOT
                figureInfoTRAJ = prepareFigure(456,size(trajs(1).Y,1),3);
            else
                figureInfoTRAJ = [];
            end
            if usejava('desktop'), figureInfoLL = prepareFigure(34,1,1); end

        end
    end
    %---------------------------------------------------------------------%
    
    
    %------------------------------- E STEP ------------------------------%
    [xpctT,XNtrpY(iEpoch),infrs] = Estep4LDS(trajs,qLDSparams,pLDSparams,...
        Nparticles,iEpoch,datadstrb,USELOGS,INFRTYPE,figureInfoTRAJ,stateinitfunc);
    % if strcmp(pLDSparams.meta, 'RandInitWithEC')
    %%%% broken: you need somehow to pull out separate C and H....
    %     xpctT = enforceNoisilyObservedControls(xpctT,C,H,pLDSparams.T);
    % end
    %---------------------------------------------------------------------%
    
    
    %------------------------ CHECK FOR CONVERGENCE ----------------------%
    fractionalImprovement = (XNtrpYold - XNtrpY(iEpoch))/abs(XNtrpYold);
    XNtrpYold = XNtrpY(iEpoch);
    fprintf('average cross entropy ')
    fprintf('= %f on EM epoch %i\n',XNtrpY(iEpoch),iEpoch);
    if usejava('desktop')
        animatePlot(figureInfoLL,2:iEpoch,XNtrpY(2:iEpoch));
    end
    if (fractionalImprovement < thr)
        fprintf('exiting with delta LL below tolerance\n');
        break
    end
    %---------------------------------------------------------------------%
    
    
    %----------------------------- M STEP --------------------------------%
    qLDSparams = Mstep4LDS(xpctT,qLDSparams,DIAGCOVS,datadstrb,USELOGS);
    %---------------------------------------------------------------------%
    
    
    
    %-------------------------- COMPARE PARAMS ---------------------------%
    similarityTransformParams(infrs,trajs,qLDSparams,pLDSparams);
    Rsqs = getNextFrameRsqs(infrs,trajs,qLDSparams,datadstrb);
    %---------------------------------------------------------------------%
    
end

qLDSparams.XNtrp = XNtrpY(iEpoch);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [xpctT,XNtrpY,infrs] = Estep4LDS(trajs,qLDSparams,pLDSparams,...
    Nparticles,iEpoch,dstrb,USELOGS,INFRTYPE,figureInfoTRAJ,stateinitfunc)
% This is a wrapper for the actual implementation of the E step.  It can
% call either the exact E-step--assuming everything is linear-Gaussian--or
% a particle-based approximation.


% BOOLEANS
BIASEDTRANSITIONS = 1;
BIASEDEMISSIONS = diff(size(qLDSparams.C)) < 0;

% some initial params (or possibly samples) may not work with the PF
ISERR = 1;
while ISERR
    try
        switch INFRTYPE
            case 'exact'
                [accumT,XNtrpY,infrs,Nsamples] = exactEstep(trajs,...
                    qLDSparams,figureInfoTRAJ);
            case 'particle'
                %%%% there will certainly be cases where pLDSparams doesn't
                %%%% have a field 'C'--because, e.g., you don't know the
                %%%% true model!
                [accumT,XNtrpY,infrs,Nsamples] = particleEstep(trajs,...
                    qLDSparams,figureInfoTRAJ,Nparticles,USELOGS,dstrb,stateinitfunc);
            otherwise
                error('not a recognized inference type -- jgm');
        end
        ISERR = 0;
    catch ME
        fprintf('Recovering from error: %s\n',ME.message);
        keyboard
        fprintf('Re-initializing samples...\n');
        if iEpoch == 1
            fprintf('...and reinitializing parameters...\n');
            qLDSparams = initLDSparams(trajs,size(qLDSparams.A,1),pLDSparams,dstrb);
        end
    end
end

% now renormalize by the number of samples in each
Ntrajs = length(trajs);
XNtrpY = XNtrpY/Ntrajs;
xpctT.mu0 = accumT.muX0/Ntrajs;
xpctT.x0x0 = accumT.X0X0/Ntrajs;
xpctT.XpXp = accumT.XpXp/(Nsamples-Ntrajs);
xpctT.XfXf = accumT.XfXf/(Nsamples-Ntrajs);
xpctT.XfXp = accumT.XfXp/(Nsamples-Ntrajs);
if BIASEDTRANSITIONS
    xpctT.XpXp = [xpctT.XpXp accumT.muXp/(Nsamples-Ntrajs);...
        accumT.muXp'/(Nsamples-Ntrajs) 1];
    xpctT.XfXp = [xpctT.XfXp accumT.muXf/(Nsamples-Ntrajs)];
end


switch dstrb
    case {'GammaFixedScale','Poisson'}
        
        % *unfortunately*, you have to pass all the data
        fprintf('not bothering to normalize sufficient stats for %s',dstrb);
        fprintf(' emissions...\n');
        if BIASEDEMISSIONS
            xpctT.Xtu  = cat(1,ones(1,Nparticles,Nsamples,...
                'like',accumT.Xtu), accumT.Xtu);
        else
            xpctT.Xtu  = accumT.Xtu;
        end
        xpctT.W = accumT.W;
        xpctT.Y = accumT.Y;
        
        
    case 'Normal'
        xpctT.XX = accumT.XX/Nsamples;
        xpctT.YX = accumT.YX/Nsamples;
        xpctT.YY = accumT.YY/Nsamples;
        if BIASEDEMISSIONS
            xpctT.XX = [xpctT.XX accumT.muX/Nsamples; accumT.muX'/Nsamples 1];
            xpctT.YX = [xpctT.YX accumT.muYX/Nsamples];
        end
    otherwise
        error('not a programmed case for E-step -- jgm');
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [accumT,XNtrpY,infrs,Nsamples] = exactEstep(trajs,...
    LDSparams,figureInfo)

% Ns
Ntrajs = length(trajs);

% init
mucov2op = @(muA,muB,Sgm)(muA*muB' + sum(Sgm,3));
if isfield(trajs,'U')
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y,...
        'controls',trajs(i).U));
else
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y));
end
if ~isfield(trajs,'SigmaYX')
    [trajs(1:Ntrajs).SigmaYX] = deal(LDSparams.SigmaYX);
end
accumT.muX0=0; accumT.muXp=0; accumT.muXf=0; accumT.muX=0; accumT.XfXp=0;
accumT.X0X0=0; accumT.XpXp=0; accumT.XfXf=0; accumT.XX=0;
accumT.muYX=0;  accumT.YY=0;  accumT.YX=0;
XNtrpY=0; Nsamples=0;

% malloc
infrs(Ntrajs,1).X = [];

% loop through cases
for iTraj = 1:Ntrajs
    
    Nsamples = Nsamples + size(trajs(iTraj).Y,2);
    
    % filter and smooth
    LDSparams.SigmaYX = trajs(iTraj).SigmaYX;
    LDSparams.T = size(trajs(iTraj).Y,2);
    filtered = KFfunc(iTraj,LDSparams);
    smoothed = RTSsmoother(LDSparams,filtered);
    %%%foo(iCase).X = RTSS.XHAT;
    
    % gather expected sufficient statistics
    % first-order
    accumT.muX0 = accumT.muX0 + smoothed.XHAT(:,1);
    accumT.muXp = accumT.muXp + sum(smoothed.XHAT(:,1:end-1),2);
    accumT.muXf = accumT.muXf + sum(smoothed.XHAT(:,2:end),2);
    accumT.muX  = accumT.muX  + sum(smoothed.XHAT,2);
    accumT.muYX  = accumT.muYX  + sum(trajs(iTraj).Y,2);
    
    % second-order
    accumT.X0X0 = accumT.X0X0 + mucov2op(smoothed.XHAT(:,1),smoothed.XHAT(:,1),smoothed.XCVRN(:,:,1));
    accumT.XpXp = accumT.XpXp + mucov2op(smoothed.XHAT(:,1:end-1),smoothed.XHAT(:,1:end-1),smoothed.XCVRN(:,:,1:end-1));
    accumT.XfXf = accumT.XfXf + mucov2op(smoothed.XHAT(:,2:end),smoothed.XHAT(:,2:end),smoothed.XCVRN(:,:,2:end));
    accumT.XfXp = accumT.XfXp + mucov2op(smoothed.XHAT(:,2:end),smoothed.XHAT(:,1:end-1),smoothed.XfXpCVRN);
    accumT.XX =   accumT.XX + mucov2op(smoothed.XHAT,smoothed.XHAT,smoothed.XCVRN);
    accumT.YX =   accumT.YX + mucov2op(trajs(iTraj).Y,smoothed.XHAT,0);
    accumT.YY =   accumT.YY + mucov2op(trajs(iTraj).Y,trajs(iTraj).Y,0);
    
    % (average) log likelihood
    XNtrpY = XNtrpY + filtered.XNtrpY;
    
    % store the whole trajectory, for plotting etc.
    infrs(iTraj,1).X = smoothed.XHAT;
    
    % plot
    if ~isempty(figureInfo)&&(iTraj==1)
        shatKF = LDSparams.C*filtered.XHATMU;
        shatRTSS = LDSparams.C*smoothed.XHAT;
        if isfield(LDSparams,'muYX')
            shatKF      = shatKF + LDSparams.muYX;
            shatRTSS    = shatRTSS + LDSparams.muYX;
        end
        y = trajs(iTraj).Y;
        
        % now write to the plot
        animatePlot(figureInfo,1:size(y,2),y,shatKF,shatRTSS);
        legend('observation','filtered','smoothed')
    end
    
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [accumT,XNtrpY,infrs,Nsamples] = particleEstep(trajs,LDSparams,...
    figureInfo,Nparticles,USELOGS,dstrb,stateinitfunc)

%%%%% TO DO
% (3) Perhaps EM can be written directly in terms of the weights, Ws, Wzz,
% rather than the expected sufficient statistics.  See Kokkala et alia.
%%%%%

% Ns
Ntrajs = length(trajs);

% initialize
z0 = stateinitfunc(Nparticles);
accumT.muX0=0; accumT.muXp=0; accumT.muXf=0; accumT.muX=0; accumT.XfXp=0;
accumT.X0X0=0; accumT.XpXp=0; accumT.XfXf=0;
switch dstrb
    case {'GammaFixedScale','Poisson'}
        accumT.Y = []; accumT.W = [];
		accumT.Xtu = zeros(size(z0,1),Nparticles,0);
    case 'Normal'
        accumT.muYX=0;  accumT.YY=0;  accumT.YX=0; accumT.XX = 0;
    otherwise
        error('not a programmed case for E-step -- jgm');
end
XNtrpY=0; Nsamples=0;

% construct the particle filter
PFfunc = setFilterFunctions(trajs,z0,LDSparams,dstrb,USELOGS);

% loop through trajectories
tic
infrs(Ntrajs,1).X = [];
for iTraj = 1:Ntrajs
    
    % filter and smooth
    [Xtu,Wf,crossEntropy] = PFfunc(iTraj);
    [Ws,thisXfXp] = particleSmoother(Xtu,Wf,LDSparams,USELOGS);
    
    
    % expected sufficient statistics for state transitions
    % second-order
    if USELOGS
        XtuW = Xtu.*shiftdim(exp(Ws),-1);
    else
        XtuW = Xtu.*shiftdim(Ws,-1);
    end
    Xbar = permute(sum(XtuW,2),[1,3,2]);
    
    accumT.X0X0 = accumT.X0X0 + XtuW(:,:,1)*Xtu(:,:,1)';
    accumT.XpXp = accumT.XpXp + sum(tensorOp(XtuW(:,:,1:end-1),...
        permute(Xtu(:,:,1:end-1),[2,1,3])),3);
    accumT.XfXf = accumT.XfXf + sum(tensorOp(XtuW(:,:,2:end),...
        permute(Xtu(:,:,2:end),[2,1,3])),3);
    accumT.XfXp = accumT.XfXp + thisXfXp;
    
    % first-order
    accumT.muX0 = accumT.muX0 + Xbar(:,1);
    accumT.muXp = accumT.muXp + sum(Xbar(:,1:end-1),2);
    accumT.muXf = accumT.muXf + sum(Xbar(:,2:end),2);
    accumT.muX  = accumT.muX  + sum(Xbar,2);
    
    
    % expected sufficient statistics for emissions
    Y = trajs(iTraj).Y;
    switch dstrb
        
        case 'GammaFixedScale'
            accumT.Y    = cat(2,accumT.Y,Y);
            accumT.Xtu  = cat(3,accumT.Xtu,Xtu);
            accumT.W    = cat(2,accumT.W,Ws);            
            %%% This is very unfortunate, but IRLS needs to recalculate the
            %%% expected sufficient statistics on each iteration (with a
            %%% new C matrix), so it needs access to *all the particles*.
            
            % for plotting
            if isfield(LDSparams,'muYX')
                outputFunc = @(XXX)(exp(LDSparams.C*XXX + LDSparams.muYX)*LDSparams.th);
            else
                outputFunc = @(XXX)(exp(LDSparams.C*XXX)*LDSparams.th);
            end
            
        case 'Poisson'
            accumT.Y    = cat(2,accumT.Y,Y);
            accumT.Xtu  = cat(3,accumT.Xtu,Xtu);
            accumT.W    = cat(2,accumT.W,Ws);            
            %%% This is very unfortunate, but IRLS needs to recalculate the
            %%% expected sufficient statistics on each iteration (with a
            %%% new C matrix), so it needs access to *all the particles*.
            
            % for plotting
            if isfield(LDSparams,'muYX')
                outputFunc = @(XXX)(exp(LDSparams.C*XXX + LDSparams.muYX));
            else
                outputFunc = @(XXX)(exp(LDSparams.C*XXX));
            end
            
            
        case 'Normal'
            % second-order
            accumT.XX = accumT.XX + sum(tensorOp(XtuW,permute(Xtu,[2,1,3])),3);
            accumT.YX = accumT.YX + Y*Xbar';
            accumT.YY = accumT.YY + Y*Y';
            
            % first-order
            accumT.muYX = accumT.muYX + sum(Y,2);
            
            % for plotting
            if isfield(LDSparams,'muYX')
                outputFunc = @(XXX)(LDSparams.C*XXX + LDSparams.muYX);
            else
                outputFunc = @(XXX)(LDSparams.C*XXX);
            end
            
            
        otherwise
            error('not a programmed case for E-step -- jgm');
    end
    
    % update the entropy
    XNtrpY = XNtrpY + crossEntropy;
    Nsamples = Nsamples + size(Xtu,3);
    
    
    % plot
    if (iTraj == 1)&&~isempty(figureInfo)
        if USELOGS
            YhatFilter   = outputFunc(permute(tensorOp(Xtu,permute(exp(Wf),[1,3,2])),[1,3,2]));
            YhatSmoother = outputFunc(permute(tensorOp(Xtu,permute(exp(Ws),[1,3,2])),[1,3,2]));
        else 
            YhatFilter   = outputFunc(permute(tensorOp(Xtu,permute(Wf,[1,3,2])),[1,3,2]));
            YhatSmoother = outputFunc(permute(tensorOp(Xtu,permute(Ws,[1,3,2])),[1,3,2]));
        end
        animatePlot(figureInfo,1:size(Y,2),trajs.Y,YhatFilter,YhatSmoother);        
        legend('obsvs','filtered','smoothed')
    end
    
    % store mean trajectories for other purposes (e.g. plotting)
    infrs(iTraj).X = Xbar;
end
toc

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function PFfunc = setFilterFunctions(obsvs,z0,qLDSparams,dstrb,USELOGS)

[Nstates,Nparticles,~] = size(z0);

% (update) the filter functions
transitionFunc = @(ZZ,t)(qLDSparams.muX + (qLDSparams.A*ZZ +...
    chol(qLDSparams.SigmaX)'*randn(Nstates,Nparticles,'like',z0)));
C = qLDSparams.C;
if isfield(qLDSparams,'muYX')
    muYX = qLDSparams.muYX;
else
    muYX = zeros(size(C,1),1);
end
switch dstrb
    case 'GammaFixedScale'
        th = qLDSparams.th;
        emissionFunc = @(YY,ZZ,t)(prod(gampdf(...
            repmat(exp(YY),[1,size(ZZ,2)]),exp(C*ZZ + muYX),th),1));
        logEmissionFunc = @(YY,ZZ,t)(logGammaFunc(YY,ZZ,C,muYX,th));
        %%% NB that the emission YY is expected to be the log of the gamma-
        %%% distributed random variable!  We use this because it's the suff
        %%% statistic with which the natural params interact linearly.
        
    case 'Poisson'
        %%%% the repmat is unfortunate and can probably be eliminated---if
        %%%% you're willing to rewrite poisspdf....
        emissionFunc = @(YY,ZZ,t)(prod(poisspdf(...
            repmat(YY,[1,Nparticles]),exp(C*ZZ + muYX)),1));
        logEmissionFunc = @(YY,ZZ,t)(sum(poisslpdf(...
            repmat(YY,[1,Nparticles]),exp(C*ZZ + muYX)),1));
    case 'Normal'
        emissionFunc = @(YY,ZZ,t)(mvnpdf(repmat(YY,[1,size(ZZ,2)])',...
            (C*ZZ + muYX)',qLDSparams.SigmaYX)');
        logEmissionFunc = @(YY,ZZ,t)(mvnlpdf(repmat(YY,[1,size(ZZ,2)])',...
            (C*ZZ + muYX)',qLDSparams.SigmaYX)');
    otherwise
        error('you never programmed this case -- jgm');
end

%%%% think about changing the dimensions of z0
if USELOGS
    PFfunc = @(ii)(particleFilter(obsvs(ii).Y,z0(:,:,ii),...
        transitionFunc,logEmissionFunc,1));
else
    PFfunc = @(ii)(particleFilter(obsvs(ii).Y,z0(:,:,ii),...
        transitionFunc,emissionFunc,0));
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function L = logGammaFunc(YY,ZZ,C,muYX,th)

% compute log{\prod_i p(y_i|z_i; C,th)}
ks = exp(C*ZZ + muYX);
logV = repmat(YY - log(th),[1,size(ZZ,2)]);
L = sum((ks-1).*logV - gammaln(ks) - exp(logV) - log(th),1);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Chat,muYhat,Ahat,muZhat] = similarityTransformParams(...
    infrs,trajs,qLDSparams,pLDSparams)
% There exists some R, muR such that X = R*Z + muR.  Hence:
%       Xnext = A*Xnow + muX + noise
%   =>  R\(Xnext - muR) = R\A*Xnow + R\muX - R\muR + R\noise
%   =>                  = R\A*(R*Znow + muR) + R\muX - R\muR + R\noise
%   =>            Znext = (R\A*R)*Znow + R\(muX + (A-I)*muR) + (R\noise)
%                       =    Ahat*Znow + muXhat              + noisehat
%       Y = C*X + muYX
%         = C*R*Z  + C*muR + muYX
%         = Chat*Z + muYhat

if isfield(trajs(1),'Z')
    
    % concatenate trajectories end to end
    Zbroad = cat(2,trajs(:).Z);
    Xbroad = cat(2,infrs(:).X);
    
    
    % regress Z (inferred) on X (true)
    [beta,Rsq,~,~,p] = linregress([Zbroad', ones(size(Zbroad,2),1)],Xbroad');
    R = beta(1:size(Zbroad,1),1:size(Xbroad,1))';
    muR = beta(size(Zbroad,1)+1,:)';
    for iCoef = 1:size(beta,2)
        fprintf('estimating parameters with R^2 = %.3f, p = %.2f\n',Rsq(iCoef),p(iCoef));
    end
    
    % state
    Ahat = R\qLDSparams.A*R;
    muZhat = R\(qLDSparams.muX + (qLDSparams.A - eye(size(qLDSparams.A)))*muR);
    
    % observations
    Chat = qLDSparams.C*R;
    muYhat = qLDSparams.C*muR;
    if isfield(qLDSparams,'muYX'), muYhat = muYhat + qLDSparams.muYX; end
    
    % say
    if isfield(pLDSparams,'A')
        fprintf('       A =              Ahat =           muX =     muXhat =\n');
        for iState = 1:max(size(Zbroad,1),size(Xbroad,1))
            fprintf('[% 0.3f, % 0.3f]   [% 0.3f, % 0.3f]    [% 0.3f]   [% 0.3f]\n',...
                pLDSparams.A(iState,1),pLDSparams.A(iState,2),...
                Ahat(iState,1),Ahat(iState,2),pLDSparams.muX(iState),muZhat(iState));
        end
        if isfield(pLDSparams,'C')
            fprintf('       C =              Chat =           muYX =     muYhat =\n');
            for iObsv = 1:min(size(Chat,1),5)
                fprintf('[% 0.3f, % 0.3f]   [% 0.3f, % 0.3f]    [% 0.3f]   [% 0.3f]\n',...
                    pLDSparams.C(iObsv,1),pLDSparams.C(iObsv,2),...
                    Chat(iObsv,1),Chat(iObsv,2),0,muYhat(iObsv))
            end
        end
    end
    
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function qLDSparams = Mstep4LDS(xpctT,qLDSparams,DIAGCOVS,dstrb,USELOGS)

USEBUILTIN = 0;


switch dstrb
    case {'GammaFixedScale','Poisson'}
        
        
        if USELOGS, W = exp(xpctT.W); else W = xpctT.W; end
                
        % use matlab's glmfit, with the option 'weights'
        if USEBUILTIN
            %%%% see if it's case sensitive for 'poisson'....
            tic
            C = qLDSparams.C;
            Y = xpctT.Y;
            for iObsv = 1:Nobsvs
                Ybig = repmat(permute(Y(iObsv,:),[1,3,2]),[1,size(Xtu,2),1]);
                C(iObsv,:) = glmfitJGM(Xtu(:,:)',Ybig(:),'poisson',...
                    'weights',W(:),'constant','off');
            end
            toc
        else 
            C = particleIRLS(W,xpctT,qLDSparams,dstrb);
        end
        
        % now update the structure
        qLDSparams = ML4LDS(xpctT,qLDSparams,'diagonal covariances',DIAGCOVS);
        if size(xpctT.Xtu,1) > size(qLDSparams.C,2)
            qLDSparams.C = C(:,2:end);
            qLDSparams.muYX = C(:,1);
        else 
            qLDSparams.C = C;
        end
        
    case 'Normal'
        
        % update parameters (M step)
        qLDSparams = ML4LDS(xpctT,qLDSparams,'diagonal covariances',DIAGCOVS);
        
        %%%% Is this necessary any longer?
        % for numerical stability
        if isa(xpctT.mu0,'gpuArray')
            Sigma2YZ = gather(qLDSparams.SigmaYX'*qLDSparams.SigmaYX);
            Sigma2X = gather(qLDSparams.SigmaX'*qLDSparams.SigmaX);
            qLDSparams.SigmaYX = gpuArray(sqrtm(Sigma2YZ));
            qLDSparams.SigmaX = gpuArray(sqrtm(Sigma2X));
        else
            qLDSparams.SigmaYX = sqrtm(qLDSparams.SigmaYX'*qLDSparams.SigmaYX);
            qLDSparams.SigmaX = sqrtm(qLDSparams.SigmaX'*qLDSparams.SigmaX);
        end
        %%%%
        
    otherwise
        error('whoops, you never programmed this case -- jgm');
        
end



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function Rsqs = getNextFrameRsqs(infrs,obsvs,qLDSparams,dstrb)

% loop through trajectories (because they may have different lengths)
SST = 0; SSE0 = 0; SSE1 = 0;
for iTraj = 1:length(infrs)
    
    % sum squared total
    Ynext = obsvs(iTraj).Y(:,2:end);
    Ycentered = Ynext - mean(Ynext,2);
    SST = SST + sum(Ycentered.^2,2);
    
    % sum squared error, R^2, using previous frame
    nextFrameError = diff(obsvs(iTraj).Y,1,2);
    SSE0 = SSE0 + sum(nextFrameError.^2,2);
    
    % sum squared error, R^2, using filter
    %%% NB that this assumes a linear backbone
    Xnext = qLDSparams.A*infrs(iTraj).X(:,1:end-1) + qLDSparams.muX;
    C = qLDSparams.C;
    if isfield(qLDSparams,'muYX'), muYX=qLDSparams.muYX; else muYX=zeros(size(C,1),1); end
    switch dstrb
        case 'GammaFixedScale'
            Yhat = (C*Xnext + muYX) + log(qLDSparams.th);
            %%% because the Y you work with is the logarithm of the gamma-
            %%% distributed variable
        case 'Poisson'
            Yhat = exp(C*Xnext + muYX);
        case 'Normal'
            Yhat = C*Xnext + muYX;
        otherwise
            error('whoops, you never programmed this case -- jgm');
    end
    nextFrameError = Ynext - Yhat;
    SSE1 = SSE1 + sum(nextFrameError.^2,2);
end
Rsqs(:,1) = 1 - SSE0./SST;
Rsqs(:,2) = 1 - SSE1./SST;


%%%% but on smoother outputs, these are actually "post-dictions"
fprintf('next-frame prediction:\n');
fprintf('Rsq prev., Rsq filter\n');
for iObsv = 1:min(size(Rsqs,1),9)
    fprintf('[% 0.3f],   [% 0.3f]\n',Rsqs(iObsv,1),Rsqs(iObsv,2));
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function z0 = sampleInitialState(qLDSparams,params,Nparticles)
% NB: (1) Assumes a normally distributed prior
%     (2) z0 has size Nstates x Nparticles x Ntrajs

Ntrajs = params.Ncases;
Ndims = params.Ndims;

if isfield(qLDSparams,'muX0')
    zmin = pinv(params.dynamics.C)*params.smin(:); % NB: for noninvertible C, this
    zmax = pinv(params.dynamics.C)*params.smax(:); %   can produce bad results!!

    x0 = sampleStatePrior(qLDSparams.muX0,qLDSparams.SigmaX0,...
        Ntrajs*Nparticles,zmin(1:Ndims),zmax(1:Ndims),0.05);
    v0 = sampleStatePrior(qLDSparams.muV0,qLDSparams.SigmaV0,...
        Ntrajs*Nparticles,[],[],0.05);
    z0 = [x0 v0];
else
    z0 = qLDSparams.mu0' + (randn(Nparticles*Ntrajs,...
        size(qLDSparams.Info0,1))*chol(qLDSparams.Info0));
    %%% mvnrnd(qLDSparams.mu0,inv(qLDSparams.Info0),Nparticles*Ntrajs);
end
z0 = reshape(z0',[size(z0,2),Nparticles,Ntrajs]);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function x = sampleStatePrior(mu,Sigma,Ncases,xmin,xmax,mrgn)

Ndims = length(mu);

if any(isinf(Sigma(:)))
    x = scalefxn(rand(Ncases,Ndims,'like',mu),...
        zeros(size(xmin),'like',mu),ones(size(xmin),'like',mu),...
        xmin+mrgn,xmax-mrgn);
elseif sum(Sigma(:)) == 0
    x = mu;
else
    x = mu' + randn(Ncases,Ndims)*chol(Sigma);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function figureInfo = prepareFigure(fignum,Ndims,Ndatavecs)

figure(fignum); clf;
hold on;
Ndims = min(Ndims,3);
for iDim = 1:Ndims
    figureInfo.subplotHandle = subplot(1,Ndims,iDim);
    hold on;
    for j = 1:Ndatavecs
        figureInfo.dataplotHandle(iDim,j) = plot(NaN,NaN);
    end
end
figureInfo.num = fignum;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function animatePlot(figureInfo,xdata,varargin)

% params
%%% Ndatavecs = length(varargin);
[Ndims,Ndatavecs] = size(figureInfo.dataplotHandle);

% pull up the current figure
set(0,'CurrentFigure',figure(figureInfo.num))
colors = 'krbcmy';
%%% try not to use any more than that
for ivec = 1:Ndatavecs
    
    ydata = gather(varargin{ivec});
    %%%Ndims = size(ydata,1);
    for iDim = 1:min(Ndims,size(ydata,1))
        set(figureInfo.dataplotHandle(iDim,ivec),'XData',xdata,'YData',...
            ydata(iDim,:),'color',colors(ivec));
    end
end
drawnow
%%% set(figureInfo.subplotHandle,'XLim',[0 Max]);
%%% if you want to keep the x or y limits fixed over the animation....

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function LDSparams = initLDSparams(trajs,Nx,datadstrb,initscheme,DIAGCOVS)
% unload
Ntrajs = length(trajs);
Y0 = arrayfun(@(iCase)(trajs(iCase).Y(:,1)),1:Ntrajs,'UniformOutput',0);
Y0 = cat(2,Y0{:});
Ybroad = cat(2,trajs(:).Y);
Ny = size(Ybroad,1);
if isfield(trajs,'U'), Ubroad = cat(2,trajs(:).U); Nu = size(Ubroad,1); end
Ts = arrayfun(@(iCase)(size(trajs(iCase).Y,2)),1:Ntrajs);
if all(Ts == Ts(1)), T = Ts(1); else T = 'variable'; end
yrclass = class(Ybroad);


switch initscheme
    
    %%% this case still only works with equal-length trajectories
    case 'withData'
        
        % initialize C using PCA (assumes no output noise)
        C = estimateOutputMatrix(Nx,Ybroad);
        LDSparams.C = C;
        
        % for building A, you'll need data that have as much rank as Nx
        xfakeshort = shortdata(Ntrajs,3,(pinv(C)*Ybroad)');
        xfakeshort = xfakeshort(:,1:rank(C),:);
        Xd = xfakeshort;
        while rank(longdata(xfakeshort)) < Nx %%% <= Nx
            Xd = diff(Xd,[],3);
            xfakeshort = cat(2,xfakeshort(:,:,1:end-1),Xd);
        end
        
        % init A based on data inferred via (initial) C
        xfakef = longdata(xfakeshort(:,:,2:end))';
        if isfield(trajs,'U')
            U = shortdata(Ntrajs,3,Ubroad);
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                longdata(U(:,:,1:end-1))';...
                ones(1,(size(xfakeshort,3)-1)*Ntrajs,yrclass)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:Nx) +...
                randn(Nx)/50;
            LDSparams.B = betaTRANS(:,(Nx+1):(Nx+Nu)) +...
                randn(Nx,Nu)/50;
        else
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                ones(1,(size(xfakeshort,3)-1)*Ntrajs,yrclass)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:end-1);
        end
        LDSparams.muX = betaTRANS(:,end);
        LDSparams.SigmaX = (xfakef*xfakef') - betaTRANS*(xfakep*xfakef');
        
        % init ICs based on initial observations
        foo = randn(Nx,Nx,yrclass)/1000;         % small & random
        LDSparams.Info0 = pinv(cov((pinv(C)*Y0)')) + foo'*foo; % so it's invertible!
        LDSparams.mu0 = mean(pinv(C)*Y0,2);
        
        % tell ML4LDS that muX may be nonzero (see NOTE below)
        %%%
        LDSparams.muX = zeros(Nx,1,yrclass);
        LDSparams.muYX = zeros(Ny,1,yrclass);
        %%%
        
        % make the output and transition noises small [wrt what??]
        foo = randn(Ny,yrclass)/10;
        LDSparams.SigmaYX = foo'*foo;
        %%% NB: this will be overwritten if LDSdata has a field SigmaYX
        
        
    case 'random'
        
        LDSparams.A = randn(Nx,Nx,yrclass);
        if isfield(trajs,'U'), LDSparams.B = randn(Nx,Nu,yrclass); end
        foo = randn(Nx,Nx,yrclass);
        LDSparams.SigmaX = foo'*foo; % 1.0e-05*[0.5 0; 0 0.5];
        LDSparams.muX = randn(Nx,1,yrclass); % [0;0];
        foo = randn(Nx,Nx,yrclass);
        LDSparams.Info0 = foo'*foo; % inv([0.0030 0; 0 0.1000]);
        LDSparams.mu0 = randn(Nx,1,yrclass); % [1.0236; 2.5];
        if Ny > Nx, LDSparams.muYX = randn(Ny,1,yrclass); end
        
        switch datadstrb
            case 'GammaFixedScale'
                %%%c = rand([1,Nx],yrclass);
                %%%LDSparams.C = repmat(c,[Nobsvs,1]);
                LDSparams.C = randn([Ny,Nx],yrclass);
                LDSparams.th = 0.05; %%% should really come in thru params
                %%%%
                LDSparams.muYX = zeros(size(LDSparams.muYX),yrclass);
                %%%%
            case 'Poisson'
                LDSparams.C = randn([Ny,Nx],yrclass);
                %%%%
                %%%LDSparams.muYX = zeros(size(LDSparams.muYX),yrclass);
                %%%%
            case 'Normal'
                LDSparams.C = randn([Ny,Nx],yrclass);
                foo = randn([Ny,Ny],yrclass)/10; %%% initially lean on obsvs
                LDSparams.SigmaYX = foo'*foo;
                %%% NB: this will be overwritten if LDSdata has a field SigmaYX
            otherwise
                error('you never programmed this case -- jgm');
        end
    
        
    case 'FactorAnalysis'
        fprintf('Initializing LDS with factor analysis\n');
        [M,bhat,~,What,SigmaXY] = FactorAnalysis(Ybroad,Nx);
        Xhat = What*(Ybroad - bhat); %  +...
        %    chol(SigmaXY)'*randn(size(M,2),size(Ybroad,2));
        %%%% 
        % The problem here is that M*What*(Y-b) + b = Y.  That is, M is the
        % pseudo-inverse of What.  Then when you fit the LDS below, the
        % emission covariance will be artificially tiny---right??
        %%%%
        
        Y = shortdata(Ntrajs,3,Ybroad');
        Z = shortdata(Ntrajs,3,Xhat');
        LDSparams = learnfullyobservedLDS(Y,Z,'diagonal covariances',DIAGCOVS);
        
        if Ntrajs < length(LDSparams.mu0)
            LDSparams.mu0 = mean(Xhat,2);
            LDSparams.Info0 = inv(cov(Xhat')');
        end
        
end

LDSparams.T = T;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function P = estimateOutputMatrix(Ncols,Y)
% Y = P*X, where size(P) = Nrows x Ncols.  Estimate P from the intrinsic
% (linear) dimensionality of Y.

% params
Nrows = size(Y,1);
NpcsMAX = 10;
%%% you could really use a better criterion for max num pcs than just 10.

% Npcs <= Nrows, Npcs <= Ncols, Npcs <= NpcsMAX
Npcs = min([NpcsMAX,Nrows,Ncols]);

% how many meaningful dimensions are there in Y?
[V,D] = eig(cov(Y'));
[~,indices] = sort(diag(D),'descend');          % sort by eigenvalue
W = V(:,indices(1:Npcs));                       % select 1st Npcs PCs
P = zeros(Nrows,Ncols,'like',Y);                % a P that does nothing
P(:,1:Npcs) = W;                                % the useful part of P

end
%-------------------------------------------------------------------------%
