function recalibration
% idea: The ostensible movement of the integrated estimated, from before to
% after recalibration, is actually the consequence of the model's visual
% space moving in the opposite direction.  So we can get an assay of the
% model's (changing) internal representation of visual space without any
% unimodal visual probe.

%%%%%%%%%%%% imaginary numbers from (presumably) Jacobians evaluated
%%%%%%%%%%%% outside their range (trig fxns).


%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -xtext -> Stest, with dimension changes
% Revised: 12/12/13
%   -changed indexing of statsL and statsN, based on changes in
%   estStatsCorePP.m.
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%



close all;

% load wts
load ../results/numhidswts/Std050

% init
stddevs = 5;
params.smpls = 15;

% generate biased testing data
[Dtest,Stest,shft] = generatebiaseddata(params,stddevs); % 10);

% adapt to the discrepant data
tic
[outBiases,stepsVWDR,wtsFinal] = integAdapt(Dtest,Stest,wts,shft,params);
toc

5

% plot the progress of the adaptation
adaptationPlot(Dtest,Stest,wtsFinal,outBiases,shft,params);

% adaptation steps
stepsEMP = diff(outBiases);

% plot empirical against theoretical, get best linear-regression fits
[beta,gamma] = plotandregress(stepsEMP,stepsVWDR,params,1:39);




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [D,S,shft] = generatebiaseddata(params,stddevs)

% get shft
shft = getshft(params,stddevs);
Nbatches = 500;

% pick a point in the middle
p0.mu = scalefxn([0.5 0.5],[0;0],[1;1],params.roboparams.thmin,params.roboparams.thmax);
p0.cov = 0;

% generate data
[D, S] = generateData(Nbatches*paramas.Ncases,params,'propbias',shft,...
    'stimulusprior',p0);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [wts,Dtrain] = train(wts,shft,params)
%%%%%%%%%
% None of this will work with the current code (02/24/16)
%   -JGM
%%%%%%%%%





% generate data
[batchdata] = generateData(200*params.Ncases,params,'propbias',shft);
batchdata = shortdata(params.Ncases,3,batchdata);

% init
iRBM = 1;
numRBMs = 1;
vishid = wts{iRBM}(1:end-1,:);
hidbiases = wts{iRBM}(end,:);
visbiases = wts{numRBMs*2-iRBM+1}(end,:)';
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
Nhid = sum(hidNums);
[Ncases,Nvis,Nbatches] = size(batchdata);
vishidinc = zeros(Nvis,Nhid);
hidbiasinc = zeros(1,Nhid);
visbiasinc = zeros(1,Nvis);

% learning parameters
epsilonw =  2e-5;
epsilonvb = 2e-5;
epsilonhb = 2e-5;
momentum = params.initialmomentum;
weightcost = params.weightcost;

% train for one epoch
for batch = 1:Nbatches,
    fprintf('.');
    
    % positive phase
    posdata = batchdata(:,:,batch);
    poshidmeans = invParamMap(posdata,vishid,hidbiases,hidDstrbs,hidNums,params);
    poshidstates = sampleT(poshidmeans,hidDstrbs,hidNums,params);
    posprods = posdata'*poshidstates;
    poshidact = sum(poshidstates);
    posvisact = sum(posdata);
    
    % negative phase
    [negvisstates, neghidstates] = CDstepper(poshidstates,vishid,...
        visbiases,hidbiases,hidDstrbs,visDstrbs,hidNums,visNums,...
        params.Ncdsteps,params);
    negprods  = negvisstates'*neghidstates;
    neghidact = sum(neghidstates);
    negvisact = sum(negvisstates);
    
    % weight/bias update
    vishidinc = momentum*vishidinc +...
        epsilonw*((posprods - negprods)/Ncases - weightcost*vishid);
    visbiasinc = momentum*visbiasinc +...
        (epsilonvb/Ncases)*(posvisact - negvisact);
    hidbiasinc = momentum*hidbiasinc +...
        (epsilonhb/Ncases)*(poshidact - neghidact);
    vishid = vishid + vishidinc;
    visbiases = visbiases + visbiasinc;
    hidbiases = hidbiases + hidbiasinc;
    
end

wts{iRBM} = [vishid; hidbiases];
wts{numRBMs*2-iRBM+1} = [vishid'; visbiases'];

[Dtrain] = longdata(batchdata);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function shatL = decoder(d,params)

% init
Nmods = length(params.mods);
Ndims = params.Ndims;
smin = params.smin;
smax = params.smax;
shatL = zeros(Ndims,Nmods);

% decode
T = displayshape(d,params);
for i = 1:Nmods
    ind = (1:Ndims) + (i-1)*Ndims;
    shatL(1,ind) = decode(T{i},[smin(:,i) smax(:,i)],params,'CoM');
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function shft = getshft(params,stddevs)

[D0, S0] = generateData(1000,params);
SINSMerr = covInCalc(D0,S0,params);
close
ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;
iDirection = 1; nDirections = 8;
phi = 2*pi*iDirection/nDirections;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
shft = sqrtm(ErrCovIn)*R*[stddevs; 0];

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [outBiases,stepsVWDR,wts] = integAdapt(Dtest,Stest,wts,shft,params)

% init
N = 40;
outBiases = zeros(N,params.Ndims*length(params.mods));
stepsVWDR = zeros(N-1,params.Ndims*length(params.mods));


outBiases(1,:) = getBiases(Dtest,Stest,wts,params);
for i=1:N-1    
    [wts,Dtrain] = train(wts,shft,params);
    stepsVWDR(i,:) = getTheorAdaptRates(Dtrain,outBiases(i,:),...
        outBiases(1,:),shft,params);
    outBiases(i+1,:) = getBiases(Dtest,Stest,wts,params);
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function outBiases = getBiases(Di,Si,wts,params)

Do = updownDBN(Di,wts,params,'Nsamples');
[statsL,statsN] = estStatsCorePP(Si,params,'CoM',Do);
outBiases =  [statsN{1}.mu(:)' statsN{2}.mu(:)']; % [b_v, b_p]

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function stepsVWDR = getTheorAdaptRates(Dtrain,outBiases,initBiases,shft,params)
% This may be overkill, but it is true that the "true" covV, expressed in
% prop space, will in general change over trials, b/c the Jacobian changes
% as vis space changes (or, equivalently, as the function mapping prop to
% vis changes).

% init
Nexamples = size(Dtrain,1);
bias = zeros(params.Ndims,length(params.mods));
bias(:,2) = shft;                               % doesn't actually get used
snglModShfts = outBiases - initBiases;
tuningCov = computetuningcovs(params);  %%%% shouldn't be in loop

% average the input covariances from the training examples
covV = 0; covP = 0;
for i = 1:Nexamples
    SILSCerr = PPCinputStats(Dtrain(i,:),tuningCov,bias);
    shatL = decoder(Dtrain(i,:),params);
    [SINSCerrMu,SINSCondCov] = SICE(shatL(:)' + snglModShfts,params,SILSCerr);
    covV = covV + squeeze(SINSCondCov(:,:,:,1));       %%% hack
    covP = covP + squeeze(SINSCondCov(:,:,:,2));       %%% hack
end
covV = covV/Nexamples;
covP = covP/Nexamples;


% theoretical adaptation rates (via VWDR)
delta = outBiases(1:2) - outBiases(3:4);        % v - p
stepsVWDR = [-delta*covV, delta*covP];


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [beta,gamma] = plotandregress(stepsEMP,stepsVWDR,params,inds)
% scatter theor. against emp. adaptation steps, w/regression lines

% init
close all
beta = zeros(2,params.Ndims*length(params.mods));
% gamma = zeros(2,params.Ndims*length(params.mods));

for i = 1:params.Ndims*length(params.mods)
    figure;
    x = stepsEMP(inds,i);
    X = [x ones(inds(end),1)];
    y = stepsVWDR(inds,i);
    beta(:,i) = (X'*X)\(X'*y);
    
    hold on;
    scatter(x,y)
    xx = linspace(min(x),max(x));
    plot(xx,xx*beta(1,i) + beta(2,i),'k')
    hold off;
    pause()
    
%     % how much is that scaling by the covariance really doing?  Try w/o it
%     figure;
%     x = stepsEMP(inds,i);
%     X = [x ones(inds(end),1)];
%     y = thing(inds,i);
%     gamma(:,i) = (X'*X)\(X'*y);
%     
%     hold on;
%     scatter(x,y)
%     xx = linspace(min(x),max(x));
%     plot(xx,xx*gamma(1,i) + gamma(2,i),'g')
%     hold off;
%     pause()
    
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function adaptationPlot(Di,Si,wts,outBiases,shft,params)

% (re)generate statistics and then plot
Do = updownDBN(Di,wts,params,'Nsamples');
[statsL,statsN] = estStatsCorePP(Si,params,'CoM',Di,Do);
ErrorStats = PPCinfo(Di,Si,statsL,statsN,params,'propbias',shft);

% call up the appropriate figure and plot the adaptations on it
figure(gcf-1)
hold on
plot(outBiases(:,1),outBiases(:,2));
hold off;
figure(gcf+1)
hold on
plot(outBiases(:,3)+shft(1),outBiases(:,4)+shft(2))
hold off

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% function shatL = getintegratedestimates(Di,wts,params)
%
% % init
% params.smpls = 15;
%
% % updown pass
% DoPRE = updownDBN(Di,wts,params,'Nsamples');
%
% % decode
% shatL = decoder(DoPRE,params);
%
% end
% %-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% function ploteffectiveshifts(Di,xi,shatOUT,shatINP,params)
% 
% nummodes = length(params.mods);
% plotinds = 1:100;
% %%%
% ntrlinds = 3:4;
% %%%
% bias = zeros(params.Ndims,nummodes);
% %%% ???
% 
% % compute tuning covariances
% tuningCov = computetuningcovs(params);
% 
% % malloc
% sopt = zeros(size(shatINP,1),size(shatINP,2)/2);
% 
% % plot
% setColors;
% figure; hold on;
% % scatter(xi(plotinds,ntrlinds(1)),xi(plotinds,ntrlinds(2)),'gx');
% % scatter(shatOUT(plotinds,ntrlinds(1)),shatOUT(plotinds,ntrlinds(2)),'go');
% for i = 1:plotinds
%     
%     % find optimal integrated estimate
%     SILSCerr = PPCinputStats(Di(i,:),tuningCov,bias);
%     [SINSCerrMu SINSCondCov] = SICE(xi(i,:),params,SILSCerr);
%     
%     
%     
%     %%% hack
%     covV = squeeze(SINSCondCov(:,:,:,1));
%     covP = squeeze(SINSCondCov(:,:,:,2));
%     covInteg = inv(inv(covV) + inv(covP));
%     
%     vhat = IK2link(shatINP(i,1:2)',params,0);
%     phat = shatINP(i,3:4)';
%     
%     %     thing = [vhat phat];
%     %     plot(thing(1,:),thing(2,:),'k');
%     
%     %     h = error_ellipse(covV,vhat,'conf',.95);
%     %     set(h,'LineWidth',1,'Color',VIScolor);
%     %     h = error_ellipse(covP,phat,'conf',.95);
%     %     set(h,'LineWidth',1,'Color',PROPcolor);
%     
%     sopt(i,:) = (covInteg*(inv(covV)*vhat + inv(covP)*phat))';
%     
%     %     h = error_ellipse(covInteg,sopt(i,:),'conf',.95);
%     %     set(h,'LineWidth',1,'Color',OPTcolor);
%     %%%
%     
%     % plot vector b/n "true" and estimated
%     thing = [xi(i,ntrlinds); shatOUT(i,ntrlinds)];
%     plot(thing(:,1),thing(:,2));
%     
%     % plot vector b/n optimal and estimated
%     thing = [sopt(i,:); shatOUT(i,ntrlinds)];
%     plot(thing(:,1),thing(:,2),'k');
%     
%     
% end
% % scatter(sopt(:,1),sopt(:,2),'k.');
% hold off;
% 
% end
% %-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% function sopt = optimalEsts(Di,xi,shatINP,shft,params)
% 
% nummodes = length(params.mods);
% bias = zeros(params.Ndims,nummodes);
% bias(:,2) = shft;
% 
% % compute tuning covariances
% tuningCov = computetuningcovs(params);
% 
% % malloc
% sopt = zeros(size(shatINP,1),size(shatINP,2)/2);
% 
% for i = 1:size(Di,1)
%     
%     % find optimal integrated estimate
%     SILSCerr = PPCinputStats(Di(i,:),tuningCov,bias);
%     [SINSCerrMu SINSCondCov] = SICE(xi(i,:),params,SILSCerr);
%     
%     
%     %%% hack
%     covV = squeeze(SINSCondCov(:,:,:,1));
%     covP = squeeze(SINSCondCov(:,:,:,2));
%     covInteg = inv(inv(covV) + inv(covP));
%     
%     vhat = IK2link(shatINP(i,1:2)',params,0) + SINSCerrMu(1,:,1)';
%     phat = shatINP(i,3:4)' + SINSCerrMu(1,:,2)';
%     
%     sopt(i,:) = (covInteg*(inv(covV)*vhat + inv(covP)*phat))';
%     %%%
%     
% end
% 
% 
% end
% %-------------------------------------------------------------------------%

