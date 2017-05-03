%%%%% This is probably mildly broken, in that the new *testing* structure
%%%%% hasn't been implemented---cf. EFH.m.  The fix should be easy.
%%%%% (02/24/16 -- jgm)


%-------------------------------------------------------------------------%
% Revised: 02/24/16
%   -altered to reflect new format of typeUnits and numsUnits
% Revised: 07/28/14
%   -for EFH trained on the full range of gains
% Created: 07/09/14
%   by JGM
%-------------------------------------------------------------------------%



clear; clc; close all;

%%
% (1) healthy
% load('results\finalwts\wtsStandard140712','wts','params');
load('results\finalwts\wts2Dallgains140725.mat','wts','params');
load('AllGainsOutputDecoder.mat','netOutputDecoder')

%%% NB!!
gains = mean([params.gmin; params.gmax]);
params.gmin = gains/4;
params.gmax = 3*gains/4;

eStatsHealthy = mastertest(wts,params,...
    'posteriorlist',{'unisensory','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);

%%
% (2) "ablated" proprioception (params.NS)
%%% params.gains = params.g*ones(1,length(params.mods));
params.gmin(strcmp('Joint-Angle',params.mods)) = 0;
params.gmax(strcmp('Joint-Angle',params.mods)) = 0;
eStatsDamaged = mastertest(wts,params,...
    'posteriorlist',{'Hand-Position','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);

%%
% (3) "recovery": train on these data (no prop) *making sure to start with
% the old weight matrix*
NepochsMax = params.NepochsMax;
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
%%% datagenargs = [datagenargs{:},{'dbndepth',iRBM,'dbnwts',wts}];
datagenargs = {'dbndepth',iRBM,'dbnwts',wts};
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(iRBM,params.numsUnits,params.typeUnits,wts,...
    params.datatype,datclass);
[rcrnts,sensory,visDstrbs] = EFHinitNs(visDstrbs,visNums,params.datatype);

% 
Nbatches = 1000;
iRBM = 1;
paramDisplay(params);


% train an EFH
params.swing = 1; %%% train on the whole range
tic; EFH; toc

% store and save the weights
wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases'];
save('results\DARPA\wtsRecovered.mat','wts','params')

% test
params.gmin = gains/4;   % to be safe
params.gmax = 3*gains/4;
eStatsRecovered = mastertest(wts,params,...
    'posteriorlist',{'Hand-Position','optimal','EFHuni'},...
    'neuralNetwork',netOutputDecoder);



%%
% (4) REPAIR
%%% train a new network, with "ICMS" on its right side
NepochsMax = params.NepochsMax;
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
%%% datagenargs = [datagenargs{:},{'dbndepth',iRBM,'dbnwts',wts}];
datagenargs = {'dbndepth',iRBM,'dbnwts',wts};
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[~,~,~,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(iRBM,params.numsUnits,params.typeUnits,wts,...
    params.datatype,dataclass);
[rcrnts,sensory,visDstrbs] = EFHinitNs(visDstrbs,visNums,params.datatype);

%
Nbatches = 1000;
iRBM = 1;
paramDisplay(params);


% The equivalent of permuting the identities of the visible units is
% permuting the weights/biases that communicate w/them
vishid = wts{1}(1:end-1,:);
hidbiases = wts{1}(end,:);
visbiases = wts{2}(end,:)';
Nmods = length(params.mods);               %%% column of the NS indices into R
Nvis = sum(params.numsUnits{1});
PROPinds = reshape(1:Nvis,[Nvis/Nmods, Nmods])*strcmp('Joint-Angle',params.mods)';
[srtd,srtinds] = sort(rand(Nvis/Nmods,1));
PROPindsPermuted = PROPinds(srtinds);
visbiases(PROPinds) = visbiases(PROPindsPermuted);
vishid(PROPinds,:) = vishid(PROPindsPermuted,:);

% restore the gains and train an EFH
params.gmin = [0  0];
params.gmax = [20 20];

datagenargs = {'dbndepth',i_rbm,'dbnwts',wts};
tic; EFH; toc

% store and save the weights
wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases'];
save('results\DARPA\wtsREPAIRed.mat','wts','params')

% test one more time
params.gmin = gains/4;      % to be safe
params.gmax = 3*gains/4; 
eStatsREPAIRed = mastertest(wts,params,...
    'posteriorlist',{'unisensory','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);


%%
% If you've cleared everything:
% Compare healthy, damaged, recovered, and repaired.
% NB that you're loading the .mat file with dates appended to their names.

%%%%
% these references to params.swing won't work with current params (1/1/17)
%%%%

clear; clc; close all;

load('wts2Dallgains140725','wts','params');
params.swing = 0.5; %%% always test on truncated range
eStatsHealthy = mastertest(wts,params,...
    'posteriorlist',{'unisensory','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);

params.gains(strcmp('Joint-Angle',params.mods)) = 0;
params.swing = 0.5; %%% always test on truncated range
eStatsDamaged = mastertest(wts,params,...
    'posteriorlist',{'Hand-Position','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);

load('results\DARPA\wtsRecovered140712.mat','wts','params')
params.swing = 0.5; %%% always test on truncated range
eStatsRecovered = mastertest(wts,params,...
    'posteriorlist',{'Hand-Position','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);

load('results\DARPA\wtsREPAIRed140712.mat','wts','params')
params.swing = 0.5; %%% always test on truncated range
eStatsREPAIRed = mastertest(wts,params,...
    'posteriorlist',{'unisensory','optimal','EFHnn'},...
    'neuralNetwork',netOutputDecoder);



%% get indices to knock out X% of the prop neurons
clear; close all; clc;
load('results\finalwts\wts2Dallgains140725.mat','wts','params');

Nmods = length(params.mods);
Npop = params.N^params.Ndims;

PROPinds = reshape(1:(Nmods*Npop),[Npop, Nmods])*...
    strcmp('Joint-Angle',params.mods)';
[srtd,srtinds] = sort(rand(Npop,1));

deadIndsOneQrtr = PROPinds(srtinds(1:(1*Npop/4)));
deadIndsOneHalf = PROPinds(srtinds(1:(2*Npop/4)));
deadIndsThreeQrtrs = PROPinds(srtinds(1:(3*Npop/4)));
deadIndsFourQrtrs = PROPinds(srtinds(1:(4*Npop/4)));

save('results\DARPA\deadinds.mat','deadIndsOneQrtr','deadIndsOneHalf',...
    'deadIndsThreeQrtrs','deadIndsFourQrtrs');




%% REPAIR for each injury
clear
%%% load('results\DARPA\wtsStandard140712','wts','params');
load('results\finalwts\wts2Dallgains140725.mat','wts','params');
load('results\DARPA\deadinds.mat','deadIndsOneQrtr','deadIndsOneHalf',...
    'deadIndsThreeQrtrs','deadIndsFourQrtrs');

deadInds = {deadIndsOneQrtr,deadIndsOneHalf,deadIndsThreeQrtrs,...
    deadIndsFourQrtrs};

NepochsMax = params.NepochsMax;
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
%%% datagenargs = [datagenargs{:},{'dbndepth',iRBM,'dbnwts',wts}];
datagenargs = {'dbndepth',iRBM,'dbnwts',wts};
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[~,~,~,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(iRBM,params.numsUnits,params.typeUnits,wts,...
    params.datatype,dataclass);
[rcrnts,sensory,visDstrbs] = EFHinitNs(visDstrbs,visNums,params.datatype);

Ncases = params.Ncases;
Nbatches = params.Nbatches;
iRBM = 1;
numRBMs = 1;



for i = 3:4
    
    paramDisplay(params);
    NepochsMax = params.NepochsMax;
    restart=0; epoch=1; erravg=0; tErravg=inf; trErravg=inf; counter=0;
    batchposhidmeans = zeros(Ncases,Nhid,Nbatches);
    setColors;
    figure(2014); clf; hold on;
    subplot(1,2,1); hold on;
    plotHandle(1) = plot(NaN,NaN);
    hold off;
    subplot(1,2,2); hold on;
    plotHandle(2) = plot(NaN,NaN);
    plotHandle(3) = plot(NaN,NaN);
    hold off;
    
    cacaca = zeros(NepochsMax,1);
    yvar = []; vvar = [];
    
    % initialize increments
    vishidinc   = zeros(Nvis,Nhid);
    hidbiasinc  = zeros(1,Nhid);
    visbiasinc  = zeros(1,Nvis);
    
    % The equivalent of permuting the identities of the visible units is
    % permuting the weights/biases that communicate w/them
    vishid = wts{1}(1:end-1,:);
    hidbiases = wts{1}(end,:);
    visbiases = wts{2}(end,:)';

    sortedDeadInds = sort(deadInds{i});
    visbiases(deadInds{i}) = visbiases(sortedDeadInds);
    vishid(deadInds{i},:) = vishid(sortedDeadInds,:);
    
    % train an EFH
    datagenargs = {'dbndepth',i_rbm,'dbnwts',wts};
    tic; EFH; toc
    
    % store and save the weights
    wts{1} = [vishid; hidbiases];
    wts{2} = [vishid'; visbiases'];
    save(['results\DARPA\wtsREPAIRed',num2str(i),'Qrtr.mat'],'wts','params')
    
    % test one more time
    eStatsREPAIRed = mastertest(wts,params);
end



%%
clear; clc;

% malloc
yrEpochs = 1:5:90;
cvrnDetEFH = zeros(length(yrEpochs),4);
cvrnDetVIS = zeros(length(yrEpochs),4);
cvrnDetICMS = zeros(length(yrEpochs),4);

for iQrtr = 1:4

    % load the neural network that was trained on the final EFH
    load(['results\DARPA\decoderREPAIRed',num2str(iQrtr),'Qrtr'],'netOutputDecoder');

    for iEpoch = 1:length(yrEpochs)
        
        % load weight/params
        wtsfile = ['wtsAtEpoch',num2str(yrEpochs(iEpoch)),'.mat'];
        load(['results\DARPA\Qrtr',num2str(iQrtr),'\',wtsfile]);
        wts = tmpwts; clear tmpwts;

        % compute three things
        params.gmin = [5  5];
        params.gmax = [15 15];
        eStats = mastertest(wts,params,'posteriorlist',{'EFHnn'},...
            'neuralNetwork',netOutputDecoder);
        cvrnDetEFH(iEpoch,iQrtr) = sqrt(det(eStats.Cvrn(:,:)));
        
        params.gmin = [5  5];
        params.gmax = [15 15];
        params.gmin(strcmp('Joint-Angle',params.mods)) = 0;
        params.gmax(strcmp('Joint-Angle',params.mods)) = 0;
        eStats = mastertest(wts,params,'posteriorlist',{'EFHnn'},...
            'neuralNetwork',netOutputDecoder);
        cvrnDetVIS(iEpoch,iQrtr) = sqrt(det(eStats.Cvrn(:,:)));
        
        params.gmin = [5  5];
        params.gmax = [15 15];
        params.gmin(strcmp('Hand-Position',params.mods)) = 0;
        params.gmax(strcmp('Hand-Position',params.mods)) = 0;
        eStats = mastertest(wts,params,'posteriorlist',{'EFHnn'},...
            'neuralNetwork',netOutputDecoder);
        cvrnDetICMS(iEpoch,iQrtr) = sqrt(det(eStats.Cvrn(:,:)));
        
    end
    
    fighandle = figure(); hold on;
    plot(yrEpochs(2:end),cvrnDetVIS(2:end,iQrtr),...
        'Color',params.VIScolor,'LineWidth',2.0)
    plot(yrEpochs(2:end),cvrnDetICMS(2:end,iQrtr),...
        'Color',params.PROPcolor,'LineWidth',2.0)
    plot(yrEpochs(2:end),cvrnDetEFH(2:end,iQrtr),...
        'Color',params.NNcolor,'LineWidth',2.0)
    xlabel('training epoch');
    hold off;
    legend off, title('');
    matlab2tikzWrapper(['errcovdet',num2str(iQrtr),'Qrtr',date],fighandle);
 
end






%% 
%%% Henid:
% If you test your trained (standard) network on data in which PROP has
% been wiped out, and (as usual) you decode in PROP space, the results are
% poor.  However, if you instead test the network on data in which *VIS*
% has been wiped out, and then (still as usual) decode in PROP space, the
% results are pretty good---very close to "optimal," which here just means
% ignoring VIS.
%
% One assumes, although it's hard to check--nonflat prior and all
% that--that this all holds if VIS and PROP are switched.  You might want
% to check this by training a network with a linear, rather than nonlinear,
% transformation b/n spaces.
%
% Now, doesn't this remind you of LMMM's or even Sam Sober's findings?
% LMMM finds biases, Sober found increased error *variance* (I think).  You
% actually find both.
%
% More explicitly: ...


%%%
% problem:
%   decoding PROP space yields subpar results when PROP is zeroed out as an
%   input, but it yields (obviously) terrible results when the model has
%   been retrained on zeroed out PROP.
% possible solutions:
%   (1) zero out VIS instead.
%       (+) easy to do conceptually (you know it will work)
%       (-) requires training a new model on VIS=0
%       (-) harder to relate to MCD results [patient goes blind??]
%       (-) actually doesn't look so bad prior to "recovery": the decoded
%       estimate is only slightly suboptimal--altho' of course optimal is
%       worse
%   (2) train a NN decoder to recover optimal estimate--either from the
%   hidden units, or from the (locally biased) center of mass on the VIS
%   population
%       (+) more congruent w/MCD to have PROP ablated
%       (-) may not work
%       (-) requires building/using NN decoder
%   [UPDATE: Works ok, but not great.  In fact, the NN doesn't improve the
%   error variance very much, altho' it does eliminate the bias.  [This is
%   somewhat surprising, since you think the (marginal) error variance is a
%   consequence of local (conditional) error variance.]
%   [UPDATE: You've got a much better NN decoder working.]
%   (3) make a datatype with uniform sampling in a rectangle of VIS space.
%       (+) you can zero out PROP and decode from VIS with no problem, so
%       congruent wih MCD work
%       (-) it'll be a pain in the ass to design such a model
%       (-) you have to train a new EFH
%   (4) try using the original weights to decode!?
