function [MedianMSEs,MSEerrorBars,xaxislabels,clrNames] = multiXprmtStats(...
    filterNames,params)
% multiXprmtStats   muliple-exeripment statistics for rEFH
%
% For some models, you will want to run multiple (identical, up to noise)
% "experiments" to produce instances, and then compute statistics across
% these experiments.

%-------------------------------------------------------------------------%
% Revised: 01/29/15
%   -more or less completely rewrote.  This is now a rational function for
%   running multiple "experiments" on filters.  Each experiment has its own
%   testing data, but the same data are used for each filter within that
%   experiment.
% Revised: 01/27/15
%   -
% Created: 01/23/15
%   by JGM
%-------------------------------------------------------------------------%

%%%%% TO DO
% (2) reinsert checks that the dynamics, params, Ncases, Nxprmts, etc., match.
% (3) relatedly, you have to set the params.machine for each loaded params
%%%%%


% compute MSEs for all models, all experiments...
eStats = multiXprmtAllStats(filterNames,params);

% ...then their medians and quartiles
Nstates = size(params.dynamics.C,2);
if any(strcmp(params.mods,'Efference-Copy'))
    Nstates = Nstates + size(params.dynamics.H,2);
end
[MedianMSEs,MSEerrorBars,xaxislabels,clrNames] =...
    acrossXprmtSummaryStats(eStats,Nstates);

 

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function eStats = multiXprmtAllStats(filterNames,params)

% init
Nxprmts = 12; %%% hard-coded
Nfilters = length(filterNames);
Ncases = params.Ncases;
T = params.dynamics.T;
yrDynamics = params.dynamics;


% which version of these models are you running?
%%% these could be made even more robust....
variant = [];
if all(size(params.dynamics.C)>=2)
    if params.dynamics.C(2,2) ~= 0
        variant = [variant,'PosVel'];
    end
end
%%% probably wrong for 2DrEFH....
if params.dynamics.A(2,1) == 0
    variant = [variant,'NoSpring'];
end


% now get all the posterior distributions
for iXprmt = 1:Nxprmts
    
    % generate fresh data for this test and training
    params.dynamics.T = T;
    params.Ncases = Ncases;
    LDSdataTest = getLDSdata(params);
    LDSdataTrain = getLDSdata(params);
    
    for iFilter = 1:Nfilters
        
        % do different things for different "models"
        filterName = filterNames{iFilter};
        switch filterName
            
            case {'EM1stOrd','EM2ndOrd','EM3rdOrd'}
                               
                filename = ['LDSparams',filterName,params.MODEL,variant,'ManyXprmts.mat'];
                load(['dynamical\finalwts\',filename]);
                if ~isequaln(rmfield(params.dynamics,{'meta'}),yrDynamics)
                    error('params and loaded file don''t match');
                else
                    name = ['EM$^',num2str(size(Allparams(iXprmt).A,2)),'$'];
                    posteriors{iFilter} = KF4PPC(LDSdataTest,Allparams(iXprmt),name);
                    clear Allparams;
                end
                
            case 'rEFH'
                filename = ['wts',params.MODEL,variant,'ManyXprmts.mat'];
                load(['dynamical\finalwts\',filename]);
                if ~isequaln(params.dynamics,yrDynamics)
                    error('params and loaded file don''t match');
                else
                    [~,~,posteriors{iFilter}] =...
                        EFHfilter(LDSdataTest,Allwts{iXprmt},params);
                end
                clear Allwts;
                
            case 'KFtrue'
                LDSparamsTrue = getLDSparams(params,'true');
                %%% a waste to do this more than once, but that's ok
                posteriors{iFilter} = KF4PPC(LDSdataTest,LDSparamsTrue,'opt');
                
            case 'sensory'
                [~,posteriors{iFilter}] = EFHfilter(LDSdataTest,[],params);
                
            case 'KFobsNoCtrl'
                [posteriors{iFilter}, ~] =...
                    obsKFwithNoInput(LDSdataTrain,LDSdataTest,params);
                
                
                %%%% currently unused
            case 'yrtest'
                % load('dynamical\nonfinalwts\wts1DrEFHwithEC150202.mat')
                % load('dynamical\finalwts\wts1DrEFHwithEC141027.mat')
                % load('dynamical\nonfinalwts\wtsHVNdamped150202.mat');
                % load('dynamical\nonfinalwts\wtsTEST150202.mat');
                % load('dynamical\nonfinalwts\NoSpringSlowest.mat','wts','params')
                load('dynamical\nonfinalwts\wts1DrEFHwithECcdfive150204.mat')
                [~,machine] = system('hostname');
                params.machine = strtrim(machine);
                % [~,~,posteriors{iFilter}] = EFHfilter(LDSdataTest,Allwts{iXprmt},params);
                [~,~,posteriors{iFilter}] = EFHfilter(LDSdataTest,wts,params);
                
            case 'KFobs'
                LDSparamsObs = getLDSparams(params,'observed',LDSdataTrain);
                posteriors{iFilter} = KF4PPC(LDSdataTest,LDSparamsObs,'obs');
                
            case 'KFobsNoCtrlDynamics'
                %%% hard-coded: 2 states, 1 state obsv., 1 ctrl, 1 ctrl obsv
                LDSparamsObsNoCtrlDynamics = learnfullyobservedLDS(LDSdataTrain,2,1,1,1);
                posteriors{iFilter} = KF4PPC(LDSdataTest,...
                    LDSparamsObsNoCtrlDynamics,'obs');
                
            otherwise
                error('no such model --- jgm');
                
        end
    end
    fprintf('...computed posteriors for experiment. %i\n',iXprmt);
    
    eStats(iXprmt,:) = testDynamics(LDSdataTest,params,0,posteriors{:});
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [MedianMSEs,MSEerrorBars,xaxislabels,clrNames] =...
    acrossXprmtSummaryStats(eStats,Nstates)


% Ns
Nmods = size(eStats,2);
Nfilters = length(eStats(1,1).N);

% malloc
MedianMSEs = zeros(Nfilters,Nmods);
MSEerrorBars = zeros(Nfilters,2,Nmods); % one positive + one negative = 2
xaxislabels = cell(Nfilters,Nmods);
clrNames = cell(Nfilters,Nmods);

for iMod = 1:Nmods
    
    
    MSEs = cell2mat(arrayfun(@(stats)(...
        (stats.Xpct(:).^2 + stats.Cvrn(:)))',...
        eStats(:,iMod),'UniformOutput',false));
    MedianMSEs(:,iMod) = median(MSEs)';
    
    % 1st and 3rd quartile, for pgfplot
    MSEerrorBars(:,:,iMod) = [prctile(MSEs,75)'-MedianMSEs(:,iMod),...
        MedianMSEs(:,iMod)-prctile(MSEs,25)'];
    
    % x-axis labels for the bar plot
    names = arrayfun(@(iFilter)(eStats(1,iMod).tags(iFilter).name),...
        1:Nfilters,'UniformOutput',false);
    %%% we can use 1 b/c this will be the same across all experiments
    xaxislabels(:,iMod) = names;
    
    % "correct" the names to something getColorName recognizes
    [names{strcmp(['EM$^',sprintf('%i',Nstates-1),'$'],names)}] = deal('EM');
    [names{strcmp(['EM$^',sprintf('%i',Nstates),'$'],names)}] = deal('EM (best)');
    
    % colors for pgfplot/LaTeX
    clrNames(:,iMod) = cellfun(@getColorName,names,'UniformOutput',false);
  
end

end
%-------------------------------------------------------------------------%