function eStats = LDStest(algorithm,Ndims,VARIANTS,SWEPTPARAM,filterNames)
% LDStest   Test filters for linear dynamical systems
%
% Writes two different tikz plots: 
%
%   'MSEsVs[SWEPTPARAM]_[MODALITY]_[DATE].tex'
%
% a simple line plot of MSEs across experiments, as for example in Fig. 2
% in LtTMSwPC; and 
%
%   'MSEsBar_[MODALITY]_[DATE].tex'
% 
% a bar plot with error bars (across experiments).  NB that the former is
% really only relevant for the SWEPTPARAM = 'Springs', 'Dampers', or
% 'Masses'; whereas the latter is only really useful SWEPTPARAM = 'Xprmts'.

%   USAGES:
%{
        filterNames = {'sensory','LDSOrd1','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{},'Xprmts',filterNames);

        filterNames = {'sensory','KFobsNoCtrl','LDSOrd2','rEFH','LDSOrd3','KFtrue'};
        LDStest('rEFH',1,{'withEC'},'Xprmts',filterNames);

        filterNames = {'sensory','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{'PosVel'},'Xprmts',filterNames);

        filterNames = {'sensory','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{'PosVel','NoSpring'},'Xprmts',filterNames);


        filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{},'Springs',filterNames);

        filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{},'Dampers',filterNames);

        filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
        LDStest('rEFH',1,{},'Masses',filterNames);

%}


%-------------------------------------------------------------------------%
% Revised: 03/04/16
%   -massively rewrote:
%       --functionized, rationalized
%       --incorporated what was formerly multiXprmtStats.m
%       --incorporated what was formerly ManyModelsResults.m
% Revised: 01/29/15
%   -pushed looping over "experiments" and model types ("filters") entirely
%   into multiXprmtStats.m.
% Revised: 01/28/15
%   -completely revised, cleaned up, to run with many experiments and error
%   bars.
% Created: 09/08/14
%   by JGM
%-------------------------------------------------------------------------%

%%%%% TO DO 
% (1) Put in checks to make sure loaded and allDynamics are the same?
%%%%%



% get the dynamics for all the different experiments
[LTIsystems,sweep,MODEL,filesuffix,params] =...
    getAllLTIsystems(algorithm,Ndims,VARIANTS,SWEPTPARAM);

% compute error statistics for all these filters
eStats = multiXprmtStats(LTIsystems,MODEL,filesuffix,filterNames,params);

% bar plot with error bars
plotAcrossXprmtSummaryStats(eStats,params);

% regular line plot of results (Figure 2 in LtTMSwPC)
plotXprmtStats(sweep,eStats,SWEPTPARAM,params.mods);



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [LTIsystems,sweep,MODEL,filesuffix,params] =...
    getAllLTIsystems(algorithm,Ndims,VARIANTS,SWEPTPARAM)
%%% NB: This assumes 1D dynamics.  It could be extend to higher dimensions
%%% with some simple conventions and indexing.

% after all, this is LDStest.m!
datatype = 'LTI-PPC';

% grab the filesuffix for loading other files
variantstr = arrayfun(@(i)(['_',VARIANTS{i}]),1:length(VARIANTS),...
    'UniformOutput',false);
MODEL = sprintf('%s_%iD_%s',algorithm,Ndims,datatype);
filesuffix = [variantstr{:},'_Many',SWEPTPARAM];

% construct the parameters
mods = {'Joint-Angle'};
SPRING = 1;
if any(strcmp(VARIANTS,'PosVel')), mods = [mods, {'Angular-Velocity'}]; end
if any(strcmp(VARIANTS,'withEC')), mods = [mods, {'Efference-Copy'}]; end
if any(strcmp(VARIANTS,'NoSpring')), SPRING = 0; end
params = setParams('datatype',datatype,'mods',mods,'SPRING',SPRING);
Nxprmts = 12;

% from the parameters file
m0 = params.dynamics.m;
dt = params.dynamics.dt;
A = params.dynamics.A;
C = params.dynamics.C;

% A is the state-transition matrix for [x_t; v_t]
k0 = -A(2,1)*m0/dt;
b0 = (1 - A(2,2))*m0/dt;

% construct the state matrices
LTIsystems.As = repmat(A,[1,1,Nxprmts]);
LTIsystems.Cs = repmat(C,[1,1,Nxprmts]);

% whatever changes across the Nxprmts experiments
switch SWEPTPARAM
    case 'Springs'
        ks = linspace(0,5,Nxprmts);
        sweep = ks;
        LTIsystems.As(2,1,:) = -ks*dt/m0;
    case 'Dampers'
        bs = linspace(0,1,Nxprmts);
        LTIsystems.As(2,2,:) = 1 - bs*dt/m0;
        sweep = bs;
    case 'Masses'
        ms = linspace(1,10,Nxprmts);
        sweep = ms;
        LTIsystems.As(2,1,:) = -k0*dt./ms;
        LTIsystems.As(2,2,:) = 1 - b0*dt./ms;
        if size(A,2) == 3, LTIsystems.As(2,3) = dt./ms; end;
    case 'Xprmts'
        sweep = 1:Nxprmts;
    otherwise
        error('unexpected case! -- jgm');
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function eStats = multiXprmtStats(LTIsystems,MODEL,filesuffix,...
    filterNames,masterparams)

% Ns
Ntraj = 40; %%% hard-coded
T = 1000;
Nxprmts = size(LTIsystems.As,3);
Nfilters = length(filterNames);
Nstates = size(masterparams.dynamics.C,2);
Ndims = masterparams.Ndims;
Nmods = length(masterparams.mods);
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end



% For each experiment...
for iXprmt = 1:Nxprmts
    
    % ...generate a set of testing data
    masterparams.dynamics.A = LTIsystems.As(:,:,iXprmt);
    masterparams.dynamics.C = LTIsystems.Cs(:,:,iXprmt);
    LDSparamsTrue = getLDSparams(masterparams.dynamics);
    
    % test data
    %%% might want to allow RTRBM, TRBM...
    if strcmp(filesuffix,'_ManyXprmts')||strcmp(filesuffix,'_withEC_ManyXprmts')
        [Rtest,Xtest,Qtest] = masterparams.getTestData(dataclass);
    else
        [Xtest,Qtest] = getLatentsLTI(Ntraj*T,T,dataclass,masterparams.NS,masterparams);
        Rtest = getDataPPC(Xtest,Qtest,masterparams);
    end
    [ShatTest,ttlSpksTest] = decodeDataPPC(Rtest,Xtest,Qtest,masterparams);
    InfoTest = GTPNposteriorInfo(ttlSpksTest,masterparams);
    unisensCmlntsTest = cumulantNeutralize(ShatTest,InfoTest,masterparams);
    multisensCmlntsTest = KFposteriorization(unisensCmlntsTest,Qtest,...
        LDSparamsTrue,masterparams);
    
    
    % %%% only create if KFobs* is among filters...........
    [Xtrain,Qtrain] = getLatentsLTI(Ntraj*T,T,dataclass,masterparams.NS,masterparams);
    Rtrain = getDataPPC(Xtrain,Qtrain,masterparams);
    ShatTrain = decodeDataPPC(Rtrain,Xtrain,Qtrain,masterparams);
    Ytrain = ShatTrain(:,:);
    
    
    % malloc
    Xpct = NaN([Ntraj*T,Ndims,Nmods,Nfilters],dataclass);
    srcs = cell(Nfilters,Nmods);
    
    % ...and for each filter...
    for iFilter = 1:Nfilters
        
        % ...get the posterior
        filterName = filterNames{iFilter};
        switch filterName
            
            case {'LDSOrd1','LDSOrd2','LDSOrd3'}
                EMfile = [filterName,'_',MODEL(end-9:end),filesuffix];
                load([getdir('data'),'RBMish/EMparams/',EMfile]);
                %%% You could enforce this check, but it would require you
                %%% to resave the Springs, Masses, Dampers files with their
                %%% respective A matrices---and maybe G matrices as well.
                % if ~isequaln(rmfield(params.dynamics,'meta'),...
                %   rmfield(masterparams.dynamics,'m'))
                %   error('LDS dynamics don''t match rEFH dynamics');
                % end
                if Nstates > size(Allparams(iXprmt).A,2)
                    [srcs{iFilter,:}] = deal('EM');
                else
                    [srcs{iFilter,:}] = deal('EM (best)');
                end
                p = KFposteriorization(unisensCmlntsTest,Qtest,...
                    Allparams(iXprmt),masterparams);
                Xpct(:,:,:,iFilter) = p.Xpct;
                clear Allparams;
                
            case 'rEFH'
                load([getdir('data'),'RBMish/EFHs/','wts_',MODEL,filesuffix]);
                masterparams.numsUnits = params.numsUnits;
                [~,~,Shat1,Info1] = testEFHPPC(Rtest,Xtest,Qtest,...
                    Allwts{iXprmt},masterparams);
                p = cumulantNeutralize(Shat1,Info1,params);
                Xpct(:,:,:,iFilter) = p.Xpct;
                [srcs{iFilter,:}] = deal('rEFH');
                clear Allwts;
                
            case 'KFtrue'
                Xpct(:,:,:,iFilter) = multisensCmlntsTest.Xpct;
                [srcs{iFilter,:}] = deal('opt');
                
            case 'sensory'
                Xpct(:,:,:,iFilter) = unisensCmlntsTest.Xpct;
                [srcs{iFilter,:}] = deal(unisensCmlntsTest.srcs{:});
                
            case 'KFobs'  
                LDSparamsObs = learnfullyobservedLDS(...
                    shortdata(Ntraj,3,Ytrain),shortdata(Ntraj,3,Xtrain));
                p = KFposteriorization(unisensCmlntsTest,Qtest,...
                    LDSparamsObs,masterparams);
                Xpct(:,:,:,iFilter) = p.Xpct;
                [srcs{iFilter,:}] = deal('obs');
                
            case 'KFobsNoCtrl'
                p = obsKFwithNoInput(Ytrain,Xtrain,unisensCmlntsTest,...
                    Qtest,masterparams);
                Xpct(:,:,:,iFilter) = p.Xpct;
                [srcs{iFilter,:}] = deal('obs');
                
            case 'KFobsNoCtrlDynamics'
                %%% hard-coded: 2 states, 1 state obsv., 1 ctrl, 1 ctrl obsv
                LDSparamsObsNoCtrlDynamics = learnfullyobservedLDS(...
                    LDSdataTrain,'alternative dimensions',[2,1,1,1]);
                p = KFposteriorization(unisensCmlntsTest,Qtest,...
                    LDSparamsObsNoCtrlDynamics,masterparams);
                Xpct(:,:,:,iFilter) = p.Xpct;
                [srcs{iFilter,:}] = deal('obs');
                %%%%%
                
            otherwise
                error('no such model --- jgm');
                
        end
    end
    
    fprintf('...computed posteriors for experiment %i\n',iXprmt);
    S = latents2stims(Xtest,Qtest.latent2stim,masterparams.mods,masterparams.Ndims);
    %%%eStats(iXprmt,:) = fltrDstrbs2ErrorStats(Xpct,S,srcs,masterparams,0);
    for iMod = 1:Nmods
        data.Xpct(1:(Ntraj*T),1:Ndims,1:Nfilters) = Xpct(:,:,iMod,:);
        data.srcs = srcs(:,iMod);
        eStats(iXprmt,iMod) = getErrorStats(data,S(:,:,iMod));
    end
    
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [MedianMSEs,MSEerrorBars,xaxislabels,clrNames] =...
    plotAcrossXprmtSummaryStats(eStats,params)

% Ns
Nmods = size(eStats,2);
Nfilters = length(eStats(1,1).N);

% malloc
MedianMSEs = zeros(Nfilters,Nmods);
MSEerrorBars = zeros(Nfilters,2,Nmods); % one positive + one negative = 2
xaxislabels = cell(Nfilters,Nmods);
clrNames = cell(Nfilters,Nmods);

for iMod = 1:Nmods
    
    % get MSEs
    MSEs = cell2mat(arrayfun(@(stats)(...
        (stats.Xpct(:).^2 + stats.Cvrn(:)))',...
        eStats(:,iMod),'UniformOutput',false));
    MedianMSEs = median(MSEs)';
    
    % 1st and 3rd quartile, for pgfplot
    MSEerrorBars = [prctile(MSEs,75)'-MedianMSEs,...
        MedianMSEs-prctile(MSEs,25)'];
    
    % x-axis labels for the bar plot
    names = arrayfun(@(iFilter)(eStats(1,iMod).tags(iFilter).name),...
        1:Nfilters,'UniformOutput',false);
    %%% we can use 1 b/c this will be the same across all experiments
    xaxislabels = names;
    
    % "correct" the names to something getColorName recognizes
    %[names{strcmp(['EM$^',sprintf('%i',Nstates-1),'$'],names)}] = deal('EM');
    %[names{strcmp(['EM$^',sprintf('%i',Nstates),'$'],names)}] = deal('EM (best)');
    
    % colors for pgfplot/LaTeX
    clrNames = cellfun(@getColorName,names,'UniformOutput',false)';
    
    % plot
    goodInds = ~isnan(MedianMSEs);
    tikzBarGraph(...
        1:sum(goodInds),...
        MedianMSEs(goodInds),...
        MSEerrorBars(goodInds,:),...
        0,...
        xaxislabels(goodInds),'estimator','mean square error','',...
        clrNames(goodInds),...
        2,0.32,'overlapping',{},...
        ['MSEsBar_',params.mods{iMod},'_',date]);
    
end



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotXprmtStats(sweep,eStats,SWEPTPARAM,mods)

for iMod = 1:length(mods)
    
    % extract the stuff for this modality, all experiments
    Xpct = cat(3,eStats(:,iMod).Xpct);
    Cvrn = cat(4,eStats(:,iMod).Cvrn);
    tags = eStats(1,iMod).tags;
    
    % get the colors and their names
    clrs = arrayfun(@(i)(getColor(tags(i).name)),...
        1:length(tags),'UniformOutput',false);
    clrs = cat(1,clrs{:});
    clrNameClrPairs = arrayfun(@(i)(...
        {getColorName(tags(i).name),getColor(tags(i).name)}),...
        1:length(tags),'UniformOutput',false);
    
    % get the MSEs...
    MSEs = arrayfun(@(iFilter)(...
        squeeze(Cvrn(:,:,iFilter,:)) + squeeze(Xpct(:,iFilter,:)).^2),...
        1:size(Cvrn,3),'UniformOutput',false);
    %%%% assumes 1D dynamics (Cvrn is a variance)
    MSEs = cat(2,MSEs{:});
    
    % ...and plot them
    figure(1545+iMod); clf;
    set(0,'DefaultAxesColorOrder',clrs)
    plot(sweep,MSEs,'LineWidth',2);
    ymin = min(MSEs(:))*0.95;
    ymax = max(MSEs(:))*1.05;    
    axis([sweep(1),sweep(end),ymin,ymax]);
    % xlabel('spring constant (kg/s^2)')
    ylabel('MSE')
    
    matlab2tikzWrapper(['MSEsVs',SWEPTPARAM,'_',mods{iMod},'_',date],...
        figure(1545+iMod),'extraColors',{clrNameClrPairs{:}});

end


end
%-------------------------------------------------------------------------%




% TO DO

% -(9) Tried recursing with V1 rather than V0---or really, V0bar.  And lo,
% it learns just as well----except that in *testing*, one must recurse with
% V0bar to get the best results (even tho' it was trained on V1!!!).
% -(12) Increase number of units??
%   [Went up to [1575 x 1350], didn't really help.  Will try even more.]
% -(13) Tried using V0 rather than V0bar in the recusion for training.
%   This seems to work just as well.
% -(14) I think you want to initialize the dRBM differently during testing.
%   Specficially, don't use the zeros vector, which is fucked up; do
%   conditional inference or something first.... [see (20)]
% (15) Use actual reach trajectories, say from JEO's data!!
% (16) Back to applying the dRBM to real data.  I think you should really
%   try using regression to decode the hidden units, rather than the
%   RBM-type method.  At least at first....
%  (21) Write a little code to *premake the training data*, and then just
%   load it in RBM rather than recreate it!
% ?(22) issue: for the EM-trained LDS, a version with an input does no
%   better than one with no input!!  And worse than one with an input
%   learned with fully observed data.
%   [see working notes (216)]
% (23) can you compute (fake) cross-entropies for the EFH?  You would need
%   to infer the cumulants of the marginal distribution of---the center of
%   mass.



% (26) re-run EM-trained model with data with higher gain (yielding
% reliabilities like those for the rEFH)
% -(27) consider replacing all Shat with Xpct (and etc.)
% (28) Consider replacing 40 trajectories of 1000 time steps apiece with 1
% trajectory of 40,000 time steps.  This could eliminate the weirdness of
% starting with a vector of all zeros/blank input for the first input.  B/c
% the minibatches are trained sequentially, you could just seed each one
% with the last one's ending vector!  ...And then you could KF on the whole
% trajectory of 40,000 time steps.
%   If you're right that the model is trying to remember all of its inputs
%   for all time, this would test that--since it presumably can't do that
%   for 40,000 time steps w/only 225 feedback units.  Of course, it can't
%   do that for 1000 time steps either, but then it isn't perfect....

% (32) see if you can explicitly store the data (Poisson, Bern) as ints??
% (33) funtionize, a la mastertest.m?
