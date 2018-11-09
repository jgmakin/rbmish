% recurrentEFHtrain
%
%   This script trains an "EFH filter" on data that encode a dynamical
%   system.  The dynamics themselves are set in setDynamics.m.  Then it
%   tests the trained network, with the classic test (error covariances
%   across all trials) as well as several visualization tools.
%
%       The *runs* are repeats over the *same* model, whereas the "experiments"
%       are repeats over (possibly) different models
%
%   You need to leave this as a script because it's often imperative to
%   break in the middle without losing data.

%-------------------------------------------------------------------------%
% Revised: 02/29/16
%   -changed to reflect new format of numsUnits, typeUnits
%   -cleaned
%   -augmented to replace both LearnDynamics.m (with which it was already
%   mostly identical) and fitLDSwithREFH.m---both of which have now been
%   retired.
% Cribbed: 05/28/15
%   -from LearnDynamics.m
%   by JGM
%-------------------------------------------------------------------------%

% init
clear; % clc;
% params = setParams('datatype','spikecounts');
params = setParams('datatype','bouncingballs');

EXPERIMENTS = 'OneOff'; % 'DifferentIntervals';
Nruns = 20;

switch EXPERIMENTS
    case 'ManyDampers'
        Nxprmts = 12;
        bs = linspace(0,1,Nxprmts);
        ks = 3*ones(1,Nxprmts);
        ms = 5*ones(1,Nxprmts);
        dt = 0.05;
        Nhids = params.numsUnits{2}*ones(1,Nxprmts);
    case 'ManySprings'
        Nxprmts = 12;
        bs = 0.25*ones(1,Nxprmts);
        ks = linspace(0,5,Nxprmts);
        ms = 5*ones(1,Nxprmts);
        dt = 0.05;
        Nhids = params.numsUnits{2}*ones(1,Nxprmts);
    case 'ManyMasses'
        Nxprmts = 12;
        bs = 0.25*ones(1,Nxprmts);
        ks = 3*ones(1,Nxprmts);
        ms = linspace(1,10,Nxprmts);
        dt = 0.05;
        Nhids = params.numsUnits{2}*ones(1,Nxprmts);
    case 'DifferentNumbersOfHids'
        %%% Nxprmts = 15;
        Nxprmts = 1;
        bs = 0.25*ones(1,Nxprmts);
        ks = 3*ones(1,Nxprmts);
        ms = 5*ones(1,Nxprmts);
        dt = 0.05;
        %%%Nsensory = params.numsUnits{1}(2);
        %%%Nhids = (1:Nxprmts)*Nsensory;
        Nhids = 900;
    case 'OneOff'
        Nxprmts = 1;
        Nruns = 1;
        Nhids = params.numsUnits{2};
        %%%%%
        % this could cause trouble in the DBN version.....
        %%%%%
    case 'DifferentIntervals'
        Nruns = 1;
        samplevec = [2,4,8,16,32,64];
        Nxprmts = length(samplevec);
        Nhids = ones(Nxprmts,1)*params.numsUnits{2};
end



% different recurrent models require slightly different params
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
EACHBATCHISATRAJ = params.EACHBATCHISATRAJ;
getLatents  = params.getLatents;
getData     = params.getData;
Ncases      = params.Ncases;
Nbatches    = params.Nbatches;
Npretrain   = params.Npretrain;
Ntest       = 5; % params.NepochsMax;
NepochsMax  = params.NepochsMax;


% Ns/malloc
Nsensory = params.numsUnits{1};
errorStatMat = zeros(Nxprmts,Nruns);
ydataTensor = zeros(params.NepochsMax/Ntest,Nxprmts,Nruns);

try
    sendmail('person@thing.edu',...
        sprintf('started training %s on %s',params.datatype,params.machine));
catch ME
    fprintf('failed to send e-mail...\n');
end

% each experiment can have *different* parameters
for iXprmt = 1:Nxprmts
    
    % for *this* experiment
    switch EXPERIMENTS
        case 'OneOff'
            % do nothing
        case 'DifferentIntervals'
            params.Nsamples = samplevec(iXprmt);
            %params = setLearningSchedules(350,200,200,'exp',params,5/params.Nsamples);
            params = setLearningSchedules(22*params.Nsamples,12*params.Nsamples,...
                12*params.Nsamples,'exp',params,110,60);
        otherwise
            params.dynamics.A = [1.0000, dt; -ks(iXprmt)/ms(iXprmt)*dt,...
                -(bs(iXprmt)/ms(iXprmt)*dt-1)];
    end
    
    
    
    params.numsUnits = {[Nhids(iXprmt),Nsensory],Nhids(iXprmt)};
    
    
    % each run has *the same* parameters
    numsUnits = params.numsUnits;
    NEFHs = length(numsUnits)-1;
    for iRun = 1:Nruns
        
        % init
        paramDisplay(params);
        wts = cell(NEFHs*2,1);
        
        % "pretraining"
        for iEFH = 1:NEFHs
            fprintf(1,'Pretraining layer %i w/EFH: %d-%d \n',...
                iEFH,sum(numsUnits{iEFH}),sum(numsUnits{iEFH+1}));
            RESTART = 1;
            
            % train
            tic; EFH; toc;
            
            % pack together weights for saving (hid => recog., vis => gener.)
            wts{iEFH} = [vishid; hidbiases];
            wts{NEFHs*2-iEFH+1} = [vishid'; visbiases'];
            
            % save, just in case
            filename = 'EncoderWtsFile';
            save(filename,'iEFH','numsUnits','wts','params','epoch');
            
        end
        ydataTensor(:,iXprmt,iRun) = gather(figmap('TestError').data);
        errorStatMat(iXprmt,iRun) = ydataTensor(end,iXprmt,iRun);
        
        
        % store these weights if they're the best
        fprintf('Finished run %i of experiment %i\n',iRun,iXprmt);
        currentStats = errorStatMat(1:iXprmt-1,:);
        currentStatVec = [currentStats(:); errorStatMat(iXprmt,1:iRun)'];
        if errorStatMat(iXprmt,iRun) == min(currentStatVec)
            bestwts = wts;
        end
        
        
        % store the weights anyway if this is a "one-off" experiment
        if strcmp(EXPERIMENTS,'OneOff'), Allwts{iRun} = wts; end
        
        save
    end
    
    
end

try
    sendmail('person@thing.edu',...
        sprintf('finished training %s on %s',params.datatype,params.machine));
catch ME
    fprintf('failed to send e-mail...\n');
end













