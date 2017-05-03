function [AllRsqsXY,AllRsqsXYX,KFRsqs,AllRsqsXYCV,AllRsqsXYXCV,KFRsqsCV] =...
    EFH4HHS(tag)
% RBM4HHS   train an EFH on HHS's data
%   This script trains a recurrent EFH on HHS's neural data corresponding
%   to a center-out reaching task.  The loaded file also includes dynamical
%   data.

%%%%%%
% This function needs to be repaired.
%%%%%% 

%-------------------------------------------------------------------------%
% Revised: 06/05/13
%   -added a line to shorten endinds (otherwise it could end up containing
%   indices of trimmed-away data)
% Revised: 06/04/13
%   -heavily reworked, cleaned up; functionized
% Created: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%


% TRAINING
% set the params that are the same for RBM1 and RBM2
params = setAllRBMparams;
params.RHtype = 'GP';
datadir = 'C:\#code\HHS\extracteddata\';

% get the training data
% load([datadir,'KFtuningdataHHS',tag]);
% [Din1,R,endinds,Xshort,Xlong,xmin,xmax,zmin,zmax,KFparams] =...
%     KFdata2XY(UnitSpikesT,St,params);
% Nstates = KFparams.Nstates;
% clear St UnitSpikesT
load(['KFdataJEO',tag]);
UnitSpikesR = UnitSpikes; clear UnitSpikes;
Sr = S; clear S;
%%%load([datadir,'KFreachingdataHHS',tag]);


[Din1,R,endinds,Zshort,Zlong,xmin,xmax,zmin,zmax,KFparams] =...
    KFdata2XY(UnitSpikesR,Sr,params);
Nstates = KFparams.Nstates;
clear Sr UnitSpikesR



% get left-half inputs...
switch params.RHtype
    case 'GB'
        % ...by training an EFH to "compress" the rates...
        [Dright,Dout1,wts1,params1,params2] =...
            trainRateEncoder(Din1,params,Nstates);
        % plotUpDown(Din1,Dout1)
    case 'GP'
        % ...or alternatively, just use the raw rates
        [Dright,wts1,params1,params2] = useRawRates(Din1,params,Nstates);
    otherwise
        error('unrecognized layer type -- jgm \n\n');
end

% train the XY associator
Din2 = cat(2,Zshort,Dright);
[Dout2,wts2,params2] = trainXYassociator(Din2,params2);
% plotUpDown(Din2,Dout2)




% TESTING
% perform conditional inference
[~,Zhatlong] = trajCondInf(Din1,Zshort,wts1,wts2,params1,params2);
ZhatUPDOWN = longdata(Dout2(:,1:KFparams.Nstates,:));
Y = R(1:size(Zlong,1),:)';
vec = (mean(Y,2)>0)';
AllRsqsXY = predictXfromY(Zlong,Y(vec,:)',ZhatUPDOWN,Zhatlong);
Xlong = scalefxn(Zlong,zmin,zmax,xmin,xmax);

% construct (four) different emission models
em = initEmissions([],...
    {Y(vec,:)',ZhatUPDOWN,Zhatlong,cat(1,Zhatlong',Y(vec,:))'});

% predict X from Y
[AllRsqsXYX,AllRsqsYX,em] = predictYfromX(Xlong,em);
KFRsqs = runFilters(Xlong,em,KFparams,endinds);
%%% figure; plot(KFRsqs'); legend('raw rates','UPDOWN','RBMed rates','both');







% CROSS VALIDATING
% load([datadir,'KFreachingdataHHS',tag]);
% [Din1,R,Xshort,Xlong,endinds] = getCVData(Sr,UnitSpikesR,...
%     xmin,xmax,zmin,zmax,KFparams,params);
% clear Sr UnitSpikesR
load(['KFtuningReachTrgDataHHS',tag]);
[Din1,R,Zshort,Zlong,endinds] = getCVData(S,UnitSpikes,...
    xmin,xmax,zmin,zmax,KFparams,params);
clear S UnitSpikes


% perform conditional inference
[~,Zhatlong] = trajCondInf(Din1,Zshort,wts1,wts2,params1,params2);
Y = R(1:size(Zlong,1),:)';
vecCV = (mean(Y,2)>0)';
AllRsqsXYCV = predictXfromY(Zlong,Y(vecCV,:)',Zhatlong);
Xlong = scalefxn(Zlong,zmin,zmax,xmin,xmax);

% construct the CV emissions
emCV = initEmissions(em([1,3,4]),...
    {Y(vec,:)',Zhatlong,cat(1,Zhatlong',Y(vec,:))'},...
    subsetting(vec,vecCV),1:Nstates,...
    [1:Nstates, Nstates + subsetting(vec,vecCV)]);

% predict X from Y
[AllRsqsXYXCV,AllRsqsYXCV,emCV] = predictYfromX(Xlong,emCV);
%%%%%%%%%%% the Rsqs are so bad that you are probably fucking up


% KF
KFRsqsCV = runFilters(Xlong,emCV,KFparams,endinds);
%%% figure; plot(KFRsqsCV'); legend('raw rates','RBMed rates','both');




%%%%%%%%%%%%
%%% (1) maybe you need to turn down the noise in the trajectories, since 
%%% they are actually "perfectly" encoded--the model seems to be treating 
%%% the data as noise......  More generally, you want to adjust the noise
%%% in sampleT.m and *perhaps* the scaling in the fxn trajectoriesToData.



% train on "Tuning," test on "Reach"
%%% raw:   0.9041    0.9122    0.7486    0.7486    0.4108    0.4146
%%% RBMed: 0.7202    0.8460    0.3817    0.5406    0.0729    0.2669
%%% both:  0.9148    0.9170    0.7594    0.7604    0.4225    0.4318 

% train on "Reach", test on (subset of) "Tuning"
%%% raw:   0.9219    0.9489    0.7788    0.8116    0.4766    0.5062
%%% RBMed  0.8934    0.9420    0.7542    0.8025    0.4712    0.4877
%%5 both:  0.9311    0.9531    0.8160    0.8397    0.5841    0.5838
%%% N = 1000


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function params = setAllRBMparams
% EFH params for all EFHs

% training params
params.NepochsMax = 90;                     % 50 for _Science_, but &c.
params.weightcost = 0.001;                  % 0.02 works well
params.counterMax = 8;
params.numtestbatches = 40;
params.Ncdsteps = 1;

% backprop params
params.BPmaxepoch = 5;
params.max_iter = 3;                        % number of linesearches
params.numMinibatches = 10;

% number of cases
params.Ncases = 20;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Din,R,endinds,Xshort,Xlong,xmin,xmax,zmin,zmax,KFparams] =...
    KFdata2XY(UnitSpikes,S,params)


% set the Kalman filter params
KFparams.Nstates = 6;
KFparams.m = 16;                              % 16 => 60 Hz (66.7 ms bins)
KFparams.dt = 1/240;                          %
KFparams = fitDynamics(S,KFparams);
KFparams.BINMETHOD = 'nonoverlappingwindow';


% bin the spike counts
[R,Xout,endinds] = binSpikeCounts(S,UnitSpikes,KFparams);

% now refit (A^16 and all that)
KFparams = refitTransitionNoise(Xout,endinds,KFparams);

% get "input" and "output" data
Din = firingRatesToData(params.Ncases,R);
[Xshort,Xlong,xmin,xmax,zmin,zmax] =...
    trajectoriesToData(Xout,params.Ncases,params.RHtype);

% trim endinds
endinds = endinds(endinds<size(Xlong,1));

	
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Dright,Dout,wts,params,params2] = trainRateEncoder(Din,params,Nstates)

params2 = params;

N = size(Din,2);
% params.numsUnits = [N, floor(N/2), floor(N/3)];
% params.numsUnits = [N, N, floor(N/3)];
params.numsUnits = {N, N};
params.typeUnits = {{'Poisson'},{'Bernoulli'}};


% training params for this EFH
params.epsilonw =  200e-5;                  % 0.02 for all three
params.epsilonvb = 200e-5;
params.epsilonhb = 200e-5;
params.initialmomentum = 0.5;
params.finalmomentum = 0.8;     % 0.9;      % 0.2 works well

% init
numsUnits = params.numsUnits;
paramDisplay(params);
numRBMs = length(numsUnits)-1;
wts = cell(numRBMs*2,1);
batchdata = Din;
[numcases,numdims,numbatches] = size(batchdata);
numvis = numdims;

% train
for i_rbm = 1:numRBMs
    fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',i_rbm,numvis,numhid);
    restart = 1;
    
    % train
    tic; EFH; toc;
    
    % pack together weights for saving (hid => recog., vis => gener.)
    wts{i_rbm} = [vishid; hidbiases];
    wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases];
    
    filename = 'EncoderWtsFile';
    save(filename,'i_rbm','numsUnits','wts','params','epoch');

    % for next time through
    numvis = numhid;
    batchdata = batchposhidmeans;    
end


% get up pass and up/down pass data
Dout = zeros(size(Din));
for iBatch = 1:numbatches
    Dout(:,:,iBatch) = updownDBN(Din(:,:,iBatch),wts,params,'means','quiet');
end



% initialize EFH number 2
Dright = batchdata;
Nhid1 = size(Dright,2);
Nhid2 = 800;
params2.N = Nstates + Nhid1;
%%% params2.numsUnits = [params2.N 100*params2.N];
params2.numsUnits = {[Nstates, Nhid1], Nhid2};
params2.typeUnits = {{'Gaussian','Bernoulli'},{'Bernoulli'}};


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotUpDown(Din,Dout)
% plot updated and unupdated data

% init
[numcases,numdims,numbatches] = size(Din);
axvec = [0 numdims min([Dout(:);Din(:)]) max([Dout(:);Din(:)])];

% plot
figure;
for iCase = 1:numcases
    iBatch = ceil(rand*numbatches);
    hold on;
    plot(squeeze(Din(iCase,:,iBatch)));
    plot(squeeze(Dout(iCase,:,iBatch)),'r');
    axis(axvec);
    hold off;
    pause;
    clf;
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Dright,wts,params,params2] = useRawRates(Din,params,Nstates)
% if using the unRBMed rates (Din), you still need to...

% ...set params2...
params2 = params;

% ...set the "first EFH" data structures to have length zero...
params.typeUnits = {};
wts = {};


params2.N = Nstates + size(Din,2);
% params2.numsUnits = [params2.N 30*params2.N];
params2.numsUnits = {[Nstates, size(Din,2)], 1000};
%%% params2.numsUnits = [params2.N 15];
%%% floor(numcases*numbatches/2/params2.N)
params2.typeUnits = {{'Gaussian','Poisson'},{'Bernoulli'}};

% let RBM2's right-half input be the unRBMed rates
Dright = Din;

% 3630 100 1000 2000 500 1500 1000

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Dout,wts,params] = trainXYassociator(batchdata,params)
% train an EFH to associate the traj's with (some version of) the rates


% training parameters
%%%%%%% These need to be fixed....
params.NepochsMax = 90;                     % 50 for _Science_, but &c.
params.epsilonw =  20e-5;
params.epsilonvb = 20e-5;
params.epsilonhb = 20e-5;
%%% lower these?? it seems to bounce around...
params.initialmomentum = 0.5;
params.finalmomentum = 0.8;     % 0.9;      % 0.2 works well
params.weightcost = 0.001;                  % 0.02 works well
params.counterMax = 8;
params.numtestbatches = 40;
params.Ncdsteps = 1;


% init
[numcases,numdims,numbatches] = size(batchdata);
iRBM = 1;
fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',i_rbm,numvis,numhid);
restart = 1;


% train
tic; EFH; toc;
    
% pack together weights for saving (hid => recog., vis => gener.)
wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases];
    
% get the up/down pass data
Dout = zeros(numcases,numdims,numbatches);
for iBatch = 1:numbatches
    Dout(:,:,iBatch) = updownDBN(batchdata(:,:,iBatch),wts,params,'means','quiet');
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Xhatshort,Xhatlong] = trajCondInf(Din,X,wts1,wts2,params1,params2)


% params
[numcases,numstates,numbatches] = size(X);

% malloc
Xhatshort = zeros(numcases,numstates,numbatches);

% run through batches...
fprintf('\n')
for i = 1:numbatches
    
    % push up to the top of the first DBN
    means = Din(:,:,i);
    for iLayer = 1:(length(params1.typeUnits)-1)
        hidDstrbs1 = params1.typeUnits{iLayer+1};
        hidNums1 = params1.numsUnits{iLayer+1};
        means = invParamMap(means,wts1{iLayer}(1:end-1,:),...
            wts1{iLayer}(end,:),hidDstrbs1,hidNums1,params1);
    end
    
    % initialize the associative memory (EFH #2)
    traj = repmat(mean(longdata(X)),numcases,1);
    brnlls = means;
    DELTA = inf;
    
    % iterate till convergence
    while norm(DELTA,'fro') > 1e-7 

        Dout = updownDBN([traj brnlls],wts2,params2,'means','quiet');
        DELTA = Dout(:,1:numstates) - traj;
        traj = Dout(:,1:numstates);
        
        % show
        if 0
            figure(132);
            imagesc(traj)
        end
    end
    
    % store
    Xhatshort(:,:,i) = traj;
    fprintf('.')
end
fprintf('\n')

% convert to long data
Xhatlong = longdata(Xhatshort);



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function AllRsqs = predictXfromY(Xlong,varargin)
% Rsq's: predict X from "Y"

% init
[M,N] = size(Xlong);
n = length(varargin);
AllRsqs = zeros(n,N);

% compute errors
for i = 1:n
    Y = varargin{i};
    [beta, AllRsqs(i,:)] = linregress([Y ones(M,1)],Xlong,'LOO');
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [AllRsqsXYX,AllRsqsYX,em] = predictYfromX(Xlong,em)
% This is a weird function: it predicts Y from X, but then measures the
% errors in predicting X, not Y, by inverting the generative distribution
% into a recognition distribution---*assuming a uniform prior over X*!

% malloc
n = length(em);
AllRsqsXYX = zeros(n,size(Xlong,2));

for i = 1:n
    [AllRsqsXYX(i,:),AllRsqsYX{i},em(i)]=fitAndInvertEmission(Xlong,em(i));
end

% Clpinv = C'*inv(SigmaYX);
% Clpinv = (C'*C)\C';
% XhatfromY = scalefxn(XhatfromY,xmin,xmax,zmin,zmax);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function D = firingRatesToData(Ncases,R)

NN = size(R,1) - mod(size(R,1),Ncases);
D = shortdata(Ncases,R(1:NN,:));

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Xnoisy,Xlong,xmin,xmax,zmin,zmax] =...
    trajectoriesToData(Xout,Ncases,typeUnits)

% current min and max
xmin = min(Xout);
xmax = max(Xout);

% scale the states into "the right range"
switch typeUnits
    case {'GB','GP'}
        zmin = xmin./sqrt(var(Xout));
        zmax = xmax./sqrt(var(Xout));
    case 'Bernoulli'
        zmin = zeros(Nstates,1);
        zmax = ones(Nstates,1);
    otherwise
        fprintf('uh-oh! --- jgm');
end
NN = size(Xout,1) - mod(size(Xout,1),Ncases);
Xlong = scalefxn(Xout(1:NN,:),xmin,xmax,zmin,zmax);

% change to (cases,dims,batches) form and possibly add noise
Xshort = shortdata(Ncases,Xlong(1:NN,:));
Xnoisy = Xshort; % + randn(size(Xtensor));


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Din1,R,Xnoisy,Xlong,endinds] = getCVData(S,UnitSpikes,...
    xmin,xmax,zmin,zmax,KFparams,params)

[R,Xout,endinds] = binSpikeCounts(S,UnitSpikes,KFparams);
Din1 = firingRatesToData(params.Ncases,R);

% change to (cases,dims,batches) form and possibly add noise
NN = size(Xout,1) - mod(size(Xout,1),params.Ncases);
Xlong = scalefxn(Xout(1:NN,:),xmin,xmax,zmin,zmax);
Xshort = shortdata(params.Ncases,Xlong);
Xnoisy = Xshort; % + randn(size(Xtensor));

% trim endinds
endinds = endinds(endinds<size(Xlong,1));

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function inds2 = subsetting(A,B)
% see subsetting.m in ../toys

inds1 = cumsum(A);                              
inds2 = unique(inds1(~A|B));

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function em = initEmissions(em,Ys,varargin)
% the varargin holds inds, if such there be.  If there aren't, then the
% other fields have not been initialized either, and conversely.

% loop through all the cells of emission data, Ys
for i = 1:length(Ys)
    em(i).Y = Ys{i};
    if isempty(varargin)
        em(i).SigmaYX = [];
        em(i).muYX = [];
        em(i).C = [];
        em(i).inds = [];
    else
        em(i).inds = varargin{i};
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Rsqs = runFilters(X,em,KFparams,endinds)


for i = 1:length(em)
    theseKFparams = KFparams;
    theseKFparams.SigmaYX = em(i).SigmaYX;
    theseKFparams.C = em(i).C;
    
    [MSE, Rsqs(:,i)] = filterAllTrials(X,em(i).Y,endinds,theseKFparams);
    
end

% to match other Rsq collections
Rsqs = Rsqs';

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [RsqXYX,RsqYX,em] = fitAndInvertEmission(X,em)
%%% maybe add similar checks for the other fields so that you can, e.g.,
%%% keep a C matrix that you filled in (by hand or whatever)


% fit
if isempty(em.SigmaYX)
    Y = em.Y;
    [beta, RsqYX, ResCV] = linregress(X,Y,'LOO');
    muYX = mean(ResCV);
    SigmaYX = cov(ResCV);
    em.inds = 1:size(Y,1);
else
    Y = em.Y(:,em.inds);
    SigmaYX = em.SigmaYX(em.inds,em.inds);
    beta = em.C(em.inds,:)';
    muYX = em.muYX(em.inds);
    Res = Y - X*beta;
    RsqYX = 1 - sum(Res.^2)./(sum(Y-repmat(mean(Y,1),size(Y,1),1)).^2);
end

% the errors should be zero mean
Y = Y - repmat(muYX,size(Y,1),1);

% get weird inversion thing
betapinv = inv(SigmaYX)*beta'*inv(beta*inv(SigmaYX)*beta');

% compute errors etc. in X space
Xhat = Y*betapinv;
Res  = Xhat - X;
MSE = sum(Res.^2)/(size(X,1)-1);
RsqXYX = 1 - MSE./var(X);

% store again
em.Y = Y;
em.SigmaYX = SigmaYX;
em.muYX = muYX;
em.C = beta';

end
%-------------------------------------------------------------------------%


