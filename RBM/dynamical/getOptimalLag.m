function [tausOpt,maxMI,Ptau,Pmax] = getOptimalLag(Z,X,lag0,lagF,params,varargin)
% getOptimalLag     Using mutual information
% USAGE:
%   [taus,maxMI,Ptau] = getOptimalLag(Z,X,params)
%   [taus,maxMI,Ptau] = getOptimalLag(Z,X,params,cntrOutInds)
%
% Find the "optimal" lag of the stimuli in X for Bernoulli probabilities X.
% That is, Z is a tensor (Ntraj x Nunits x T) of conditional probabilities,
% p(Z=1|r), for Bernoulli random variables, and X is a tensor (Ntraj x T)
% of the corresponding stimuli.  This function sample-averages out R under
% p(r|s) to get p(Z=1|s), but it only reports the result, Ptau, at the lag,
% taus, that maximizes the mutual information between X and Z.  This MI is
% also reported as maxMI.
%
% The usual params structure is also required [[but only for something
% that's current unused!]]

%-------------------------------------------------------------------------%
% Revised: 03/24/15
%   -changed loop order, relata: massively sped up
%   -added shuffling to control for significance of MIs
%   -made lag0 and lagF arguments to the function, rather than fixed params
% Revised: 03/20/15
%   -added varargin for the case of trajectories of non-constant size (i.e.
%   length in samples).  See extractCenterOutInds.m for details on this
%   argument.
% Revised: 03/18/15
%   -functionized main loop contents
% Revised: 03/09/15
%   -massively revised:
%   -allowed future lags
%   -now also returns Pmax = max_x[  P(Z_t=1|X_{t-tau})  ], the maximum of
%   the likelihood (which is not the MLE!)
% Created: 07/20/14
%   by JGM
%-------------------------------------------------------------------------%


% params
Nunits = size(Z,2);
Nxbins = 30;
Nreshuffles = 20;
pvalue = 0.95; % 0.99;


% get edge vector
edgevector = getEdgeVector(X,Nxbins,params,'3STD');
lagvector = lag0:lagF;
Nlags = length(lagvector);

% malloc (order: Nunits x Nlags x Nxbins)
P = zeros(Nxbins,Nunits,Nlags);     % P(Z_t=1|X_{t-tau}), all tau, units
MI = zeros(Nunits,Nlags);           % MI(X_{t-tau},Z_t), all tau, all units
MIshuffled = zeros(Nreshuffles,Nunits);
thr = zeros(Nunits,Nlags);

% get MI(Z;X) and P(Z=1|X)
tic;
fprintf('looping through lags...\n');
for iLag = 1:Nlags
    
    % ...
    [z,x] = siftData(Z,X,lagvector(iLag),varargin{:});
    
    % Pr(Z=1|X), I(Z;X)
    [P(:,:,iLag), MI(:,iLag)] = getPrAis1givenB(z,x(:),edgevector);
    
    % I(Zshuffled;X), shuffle over time, to get significance thresholds
    for iShuffle = 1:Nreshuffles
        [~,randomizedInds] = sort(rand(size(z)));
        randomizedInds = randomizedInds + ((1:Nunits)-1)*size(z,1);
        zShuffled = z(randomizedInds);
        [~,MIshuffled(iShuffle,:)] = getPrAis1givenB(zShuffled,x(:),edgevector);
    end
    thr(:,iLag) = quantile(MIshuffled,pvalue);
    
    fprintf('.');
end
fprintf('\n');
toc

% get the MI-maximizing lag, and the MI(Z;X) and P(Z=1|X) at these lags
MI = MI.*(MI>thr);              % kill MI below the threshold
[maxMI,indMax] = max(MI,[],2);
tausOpt = lagvector(indMax);

Ptau = arrayfun(@(iUnit)(P(:,iUnit,indMax(iUnit))),1:Nunits,'UniformOutput',false);
Ptau = cat(2,Ptau{:})';
Pmax = squeeze(max(P));         % max_x [P(Z_t=1|X_{t-tau})], all tau,units

%  figure(7); clf; imagesc(P'); pause() %%% cf. Fig. 2C in Mulliken2008
[~,bbb] = sort(tausOpt);
figure(5346); clf; 
imagesc(MI(bbb,:)./max(MI(bbb,:),[],2)); 
axis xy
%%% cf. Fig. 3A in Mulliken2008.

% find lag that maximizes mutual information
% [maxMI,indMax] = max(MI,[],2);
% taus = indMax-1;


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [zSifted, xSifted] = siftData(Z,X,lag,varargin)

% Ns
[Ntraj,Nunits,T] = size(Z);

if isempty(varargin)
    % indices
    respTimes = intersect((1:T)-lag, 1:T);
    stimTimes = intersect((1:T)+lag, 1:T);
    
    zSifted = longdata(Z(:,:,respTimes));
    xSifted = X(:,stimTimes)';
    xSifted = xSifted(:);
else
    cntrOutInds = varargin{1};
    
    respTimes = [];
    stimTimes = [];
    for iTraj = 1:length(cntrOutInds)
        [rt,st] = arrayfun(@(iInd)(...
            siftInds(cntrOutInds{iTraj}(:,iInd),lag,T,Ntraj,iTraj)),...
            1:size(cntrOutInds{iTraj},2),'UniformOutput',false);
        respTimes = [respTimes,rt{:}];
        stimTimes = [stimTimes,st{:}];
    end
    
    % store
    zSifted = zeros(length(respTimes),Nunits);
    for iUnit = 1:Nunits
        thisZ = squeeze(Z(:,iUnit,:))';
        zSifted(:,iUnit) = thisZ(respTimes);
    end
    X = X';
    xSifted = X(stimTimes);    
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [respTimes, stimTimes] = siftInds(cntrOutInds,lag,T,Ntraj,iTraj)

a = sub2ind([T,Ntraj],cntrOutInds(1),iTraj);
b = sub2ind([T,Ntraj],cntrOutInds(2),iTraj);
if lag>0
    respTimes = a:(b-lag);
    stimTimes = (a+lag):b;
else
    respTimes = (a-lag):b;
    stimTimes = a:(b+lag);
end

end
%-------------------------------------------------------------------------%
