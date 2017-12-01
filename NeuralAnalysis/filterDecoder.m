function [XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,Bv,LDSparams] =...
    filterDecoder(Vtrain,Xtrain,Vtest,Xtest,Ntraj,SStot,Nmsperbin)
% filterDecoder    Decode variables from a latent-state filter
% 
% USAGE:
%   [XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,Nv,LDSparams] =...
%       filterDecoder(Vtrain,Xtrain,Vtest,Xtest,Ntraj,SStot);
% 
% Suppose you fit to observed dynamical data a latent-variable model--say, 
% an LTI system, or an rEFH--but you also have in hand some other set of 
% latent variables, X, that you'd like to "decode." You could simply decode
% X from the variable V of the latent-variable model (its hidden state or
% perhaps even more).  But you could furthermore use those decoded
% estimates as themselves observations from an LTI system with X on the
% backbone--and therefore apply a Kalman filter.  This function acquires
% both of these decoders, static (Bv) and dynamic (LDSparams), on training 
% pairs Vtrain and Xtrain, and returns estimates (XhatStatic,XhatDynamic) 
% and coefficients of determination (RsqStatic,RsqDynamic) for the testing
% pairs, Vtest,Xtest. 
%
% The number of "trajectories" (Ntraj) which compose the test data, and
% their total sum of squares (SStot), must also be passed as inputs.

%-------------------------------------------------------------------------%
% Revised: 03/13/17
%   -cleaned up, generalized so that it applies to LDS models as well as
%   rEFHs.
%   -renamed: linearREFHDecoder ->
% Revised: 03/11/17
%   -finally settled on a regularization/decoding scheme
% Revised: 01/11/17
%   -part of the Grand Revision
% Revised: 09/27/16
%   -added KF decoder based on observed LDS model
%   -functionized, rationalized, cleaned up
% Created: ??/??/16
%   by JGM
%-------------------------------------------------------------------------%

%%%%%% bring this in from the outside (when you're ready)
USELAGS = 0;
USEUKF = 0;
%%%%%%

% train decoders, apply static, apply dynamic

if USEUKF
    [Bv,LDSparams,lags] = getLatentVariableDecoders_Li(Vtrain,Xtrain,Ntraj,Nmsperbin);
    [XhatStatic,RsqStatic] = testStaticDecoder(Vtest,Xtest,Bv,lags,SStot);
    [XhatDynamic,RsqDynamic] = testDynamicDecoder_Li(Vtest,Xtest,LDSparams,SStot);
else
    [Bv,LDSparams,lags] = getLatentVariableDecoders(Vtrain,Xtrain,Ntraj,USELAGS);
    [XhatStatic,RsqStatic] = testStaticDecoder(Vtest,Xtest,Bv,lags,SStot);
    [XhatDynamic,RsqDynamic] = testDynamicDecoder(XhatStatic,Xtest,LDSparams,SStot);
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Bvx,LDSparamsObs,lags] = getLatentVariableDecoders(V,X,Ntraj,USELAGS)

% Ns
[Nsamples,Nv] = size(V);
Ndims = 2;                                      % **assume** a plane
Nstates = size(X,2);
Norders = Nstates/Ndims;

% malloc/init
lags = zeros(Norders,1);

% fit static decoder
if USELAGS
    Bvx = zeros(Nv,Nstates,'like',X);
    Xhat = zeros(Nsamples,Nstates,'like',X);
    for iOrder = 1:Norders
        inds = (iOrder-1)*Ndims + (1:Ndims);
        [Bvx(:,inds),XhatBest,lags(iOrder)] = getLaggedDecoder(V,X(:,inds));
        %%% you don't actually need to pad here, because you just throw it
        %%% out below....
        Xhat(:,inds) = [XhatBest(1:lags(iOrder),:); XhatBest];
        fprintf('using %i lags for derivative %i...\n',...
            lags(iOrder),iOrder-1);
    end
    
    % don't bother to learn on the pre-padded observations
    Xhat = Xhat((max(lags)+1):end,:);
    X    = X((max(lags)+1):end,:);
    
else
    [Bvx,~,~,Xhat] = linregress(V,X,'L2 regularize',(Nv/Nsamples)^(10/3));
end

% fit dynamic decoder
LDSparamsObs = learnfullyobservedLDS(shortdata(Ntraj,3,Xhat),...
    shortdata(Ntraj,3,X));
LDSparamsObs.mu0 = mean(X)';                    % use *all* data to get 
LDSparamsObs.Info0 = inv(cov(X));               %   initial state


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Bvx,LDSparamsObs,lags] = getLatentVariableDecoders_Li(V,X,Ntraj,Nmsperbin)

% Ns
[Nsamples,Nv] = size(V);
ptaps = 1;%%%floor(128/Nmsperbin) + 1;
ftaps = 1;%%%floor(128/Nmsperbin);
Ndims = 2;                                      % **assume** a plane
Nstates = size(X,2);
Norders = Nstates/Ndims;
lags = zeros(Norders,1);

% static decoder
Bvx = linregress(V,X,'L2 regularize',(Nv/Nsamples)^(10/3));

% dynamic decoder (UKF)
LDSparamsObs = learnfullyobservedLDS_Li(shortdata(Ntraj,3,V),...
    shortdata(Ntraj,3,X),'ftaps',ftaps,'ptaps',ptaps);
LDSparamsObs.mu0 = mean(X)';                    % use *all* data to get
LDSparamsObs.Info0 = inv(cov(X));               %   initial state

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [XhatStatic,RsqStatic] = testStaticDecoder(V,X,Bv,lags,SStot)

% decode
XhatStatic = V*Bv;

% apply lags
if any(lags>0)
    Ndims = 2;                                  % *assume* a plane
    Nstates = size(X,2);
    Norders = Nstates/Ndims;
    for iOrder = 1:Norders
        %%% This is a little unforunate, because the positions, velocities,
        %%% and accelerations that you're padding may be mutually
        %%% inconsistent
        inds = (iOrder-1)*Ndims + (1:Ndims);
        XhatStatic(:,inds) = [...
            XhatStatic(1:lags(iOrder),inds);...
            XhatStatic(1:(end-lags(iOrder)),inds)];
    end
end

% calculate coefficients of determination
RsqStatic  = 1 - sum((XhatStatic - X).^2)./SStot;
fprintf('R^2 for static map, [composite] -> X:\n');
rprintr(RsqStatic); fprintf('\n');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [XhatDynamic,RsqDynamic] = testDynamicDecoder(Xhat,Xtest,...
    LDSparams,SStot)

% test Kalman filter on LDS
LDSparams.T = size(Xhat,1);
KFdstrbs = KalmanFilter(LDSparams,Xhat');
XhatDynamic = KFdstrbs.XHATMU';
%%% This would be cheating, of course, but it could be interesting:
% RTSSdstrbs = RTSsmoother(LDSparamsObs,KFdstrbsTrainObs);
% XhatZhat = RTSSdstrbs.XHAT';
%%%

% calculate coefficients of determination
RsqDynamic = 1 - sum((XhatDynamic - Xtest).^2)./SStot;
fprintf('R^2 for dynamic (KF) map, [decoded X] -> X:\n');
rprintr(RsqDynamic); fprintf('\n');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function rprintr(Rsqs)

arrayfun(@(ii)(fprintf('% 0.2f ',Rsqs(ii))),1:length(Rsqs));

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function UKFparams = learnfullyobservedLDS_Li(Vtrain,Xtrain,varargin)

% taps, etc.
fparams = [];  % extra params for f, not needed
ftaps = defaulter('ftaps',2,varargin{:});       % future taps
ptaps = defaulter('ptaps',3,varargin{:});       % past taps
mtaps = defaulter('mtaps',1,varargin{:});       % backbone; <=ftaps+ptas
htaps = defaulter('htaps',0,varargin{:});       % spiking history taps
alphamax = 16; % maximum regularization (ridge regression) parameter

% f: backbone (px py vx vy ax ay) -> variables encoded linearly by spikes
f = @(a,b)(a);
% f = @(a,b)([...
%     a;...
%     sqrt(a(1,:).^2+a(2,:).^2);...
%     sqrt(a(3,:).^2+a(4,:).^2);...
%     sqrt(a(5,:).^2+a(6,:).^2);...
%     a(1,:).*a(3,:);...
%     a(2,:).*a(4,:);...
%     ]);
%%%

% hack
[Xtrain,Vtrain] = longdata(gather(Xtrain),gather(Vtrain));
%%% Li's UKF doesn't accept gpuData or multiple trajectories. Luckily, you 
%%% know there is only one.

% fit
UKFparams = fit_ar_ukf_hist_rrcv(Xtrain,Vtrain,f,fparams,...
    ftaps,ptaps,mtaps,htaps,alphamax);

% re-fit initial states from *all* training data
UKFparams.xinit = mean(Xtrain);
UKFparams.vinit = cov(Xtrain);
UKFparams.ptaps = ptaps;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [XhatDynamic,RsqDynamic] = testDynamicDecoder_Li(Vtest,Xtest,...
    UKFparams,SStot)

% test
XhatDynamic = ar_ukf_hist(gather(Vtest(1:100,:)),UKFparams,...
    gather(UKFparams.xinit),gather(UKFparams.vinit),1);


% %%%%%%
LDSparams.A = UKFparams.F;
LDSparams.muX = zeros(length(UKFparams.xmean)*UKFparams.taps,1);
LDSparams.SigmaX = UKFparams.Q;
LDSparams.C = UKFparams.H.*UKFparams.xastd';        % see lines 163-164 of 
LDSparams.muYX = -LDSparams.C*UKFparams.xamean;     %   ar_ukf_hist.m
LDSparams.SigmaYX = UKFparams.R;
xinit = (UKFparams.xinit(:) - UKFparams.xmean)./UKFparams.xstd;
vinit = (UKFparams.vinit./UKFparams.xstd)./UKFparams.xstd';
LDSparams.mu0 = repmat(xinit,UKFparams.taps,1);
LDSparams.Info0 = kron(eye(UKFparams.taps),inv(vinit));


% now standardize the observations (they were fit that way)
temp = (Vtest - UKFparams.ymeans')./UKFparams.ystd';
LDSparams.T = size(temp,1);

% run (about 30% faster off the GPU)
tic
KFdstrbs = KalmanFilter(LDSparams,gather(temp(1:100,:))','lightweight',1);
XhatDynamic2 = (KFdstrbs.XHATMU.*repmat(UKFparams.xstd,UKFparams.taps,1) +...
    repmat(UKFparams.xmean,UKFparams.taps,1))';
toc





% put back in the right shape
Nstates = size(Xtest,2);
indsNow = (1:Nstates) + (UKFparams.ptaps-1)*Nstates;
XhatDynamic2 = XhatDynamic2(:,indsNow);
% %%%%%


% compute R^2 using predictions for *current* kinematic vars
Nstates = size(Xtest,2);
indsNow = (1:Nstates) + (UKFparams.ptaps-1)*Nstates;
XhatDynamic = XhatDynamic(:,indsNow);
RsqDynamic = 1 - sum((XhatDynamic - Xtest).^2)./SStot;


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Bbest, XhatBest, lag] = getLaggedDecoder(V,X)

[Nsamples,Nv] = size(V);
RsqsOld = [-inf,-inf];
lag = 0;
while true
    [Btry,Rsqs,~,Xhat] = linregress(V(1:(end-lag),:),X((lag+1):end,:),...
        'L2 regularize',(Nv/Nsamples)^(10/3));
    if any(Rsqs < RsqsOld), lag = lag-1; break; end
    Bbest = Btry;
    XhatBest = Xhat;
    RsqsOld = Rsqs;
    lag = lag+1;
end

end
%-------------------------------------------------------------------------%