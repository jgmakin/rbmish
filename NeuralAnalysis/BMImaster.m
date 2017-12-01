% BMImaster
% Loop through all JEO and Kording BMI data (as of now), and train decoders

%-------------------------------------------------------------------------%
% Revised: 03/13/17
%   -massively.  The training and testing guts are now in
%   filtersForNeuralData.m.  This script can loop across bin sizes as well
%   as sessions.
% Created: 01/21/17
%   by JGM
%-------------------------------------------------------------------------%

% init
clc; clear;
if checkGPUavailability
    dataclass = 'gpuArray'; 
else
    dataclass = 'double';
end
Nreachesmin = 300;

%--------------------- PARAMS TO BE SET BY THE USER ----------------------%
decodernames = {'static','kfobs','kfemstatic','kfemdynamic','ukf',...
    'refhstatic','refhdynamic','ukfwithtaps','wf'};
sensory_unit_type = 'Poisson'; 
swept_param = 'hidsensoryratios';        % see sweepable_params below
monkeys = {'Indy','Loco','Chewie'};
%%% You can't currently train the StandardNormal and Poisson models in the
%%% same run of this script.
TRAINNEWREFHS = 1;
TRAINNEWLDS = 0;
Nstates = 6;
TOPLOT = 0;
%-------------------------------------------------------------------------%



% set defaults for swept and non-swept parameters
sweepable_params =...
    {'binwidths','trainingtimes','fracneurons','hidsensoryratios','phidtargets'};
switch swept_param
    case 'binwidths'
        binwidths = [16,32,64,128];
        trainingtimes = 320;
        fracneurons = 1;
        hidsensoryratios = 4;
        phidtargets = 0.05;
        Nswept = length(binwidths);
    case 'trainingtimes'
        binwidths = 64;
        trainingtimes = [80, 160, 320];
        fracneurons = 1;
        hidsensoryratios = 4;
        phidtargets = 0.05;
        Nswept = length(trainingtimes);
    case 'fracneurons'
        binwidths = 64;
        trainingtimes = 320;
        fracneurons = (1/2).^(3:-1:0);
        hidsensoryratios = 4;
        phidtargets = 0.05;
        Nswept = length(fracneurons);
    case 'hidsensoryratios'
        binwidths = 64;
        trainingtimes = 320;
        fracneurons = 1;
        hidsensoryratios = [1/2,1,2,3,4];
        phidtargets = 0.05;
        Nswept = length(hidsensoryratios);
    case 'phidtargets'
        binwidths = 64;
        trainingtimes = 320;
        fracneurons = 1;
        hidsensoryratios = 4;
        phidtargets = [0.025, 0.05, 0.1, 0.2, 0.4];
        Nswept = length(phidtargets);
    otherwise
        error('unrecognized parameter for sweeping -- jgm')
end
sweep_incrementer = @(iS)(...
    strcmp(sweepable_params, swept_param)*(iS-1) +...
    ones(size(sweepable_params)));
sweep_name = [upper(swept_param(1)), swept_param(2:end-1), 'Sweep'];

% vector lengths
Nmonkeys = length(monkeys);
Ndecoders = length(decodernames);


% loop over monkeys
for iMonkey = 1:Nmonkeys
    monkey = monkeys{iMonkey};

    localdir = [monkey,'_datafiles',filesep];
    clientloc = [getdir('data'),localdir];
    spreadsheet = [monkey,'_HC_Sessions.csv'];

    % unfortunately, you have to know these
    fieldnames = {'date','seqnum','# reaches','area(s)','target grid','minnie?'};
    fileID = fopen([getdir('code'),'NeuralAnalysis',filesep,spreadsheet],'r+');
    C = textscan(fileID,'%s %f %f %s %s %c','Delimiter',',');
    fclose(fileID);
    
    % on minnie and at least Nreachesmin reaches
    goodSessions = (C{strcmp(fieldnames,'minnie?')} == 'y')&...
        (C{strcmp(fieldnames,'# reaches')} >= Nreachesmin);
    
    % the relevant info for loading, just for the "good" sessions
    dates   = C{strcmp(fieldnames,'date')}(goodSessions);
    seqnums = C{strcmp(fieldnames,'seqnum')}(goodSessions);
    mods    = C{strcmp(fieldnames,'area(s)')}(goodSessions);
    
    % malloc
    Nsessions = sum(goodSessions);
    Rsqs = zeros(Nsessions,Nstates,Nswept,Ndecoders,dataclass);
    NdataTest = zeros(Nsessions,1,Nswept,dataclass);
    
    
    % loop across sessions and bins
    for iSession = 1:Nsessions
        
        % load or download data
        session = dates{iSession};
        year    = session(1:4);
        month   = session(5:6);
        day     = session(7:8);
        switch monkey
            case 'Jackson'
                filename = sprintf('spikes_and_kinematics_%s_%02d.mat',...
                    session,seqnums(iSession));
            case 'Chewie'
                filename = sprintf('%s_%s%s%s.mat',monkey,month,day,year);
            otherwise
                filename = sprintf('spikes_and_kinematics_%s_%02d.mat',...
                    session,seqnums(iSession));
        end
        %serverloc = sprintf('../../sabes_data2/joeyo/data/%s/%s/%s/%s/processed/',...
        %    lower(monkey),year,month,session);
        serverloc = sprintf('data/%s_datafiles/',monkey);
        getRemoteData(filename,serverloc,clientloc);
        filesuffix = sprintf('%s%s%s_%02d.mat',year,month,day,seqnums(iSession));
        
        % eliminate initial positions outside of target space
        if 0
            load([clientloc,filename],...
                'chan_names','cursor_pos','finger_pos','spikes','t','target_pos');
            i0 = find(prod([cursor_pos > min(target_pos),...
                cursor_pos < max(target_pos)],2),1);
            t = t(i0:end);
            cursor_pos = cursor_pos(i0:end,:);
            finger_pos = finger_pos(i0:end,:);
            target_pos = target_pos(i0:end,:);
            save([clientloc,filename],...
                'chan_names','cursor_pos','finger_pos','spikes','t','target_pos');
        end
        
        
        % sweep the parameter
        for iSweptparam = 1:Nswept
            sweepable_params_inds = sweep_incrementer(iSweptparam);
            this_binwidth = binwidths(sweepable_params_inds(...
                strcmp(sweepable_params,'binwidths')));
            this_trainingtime = trainingtimes(sweepable_params_inds(...
                strcmp(sweepable_params,'trainingtimes')));
            this_fracneurons = fracneurons(sweepable_params_inds(...
                strcmp(sweepable_params,'fracneurons')));
            this_hidsensoryratio = hidsensoryratios(sweepable_params_inds(...
                strcmp(sweepable_params,'hidsensoryratios')));
            this_phidtarget = phidtargets(sweepable_params_inds(...
                strcmp(sweepable_params,'phidtargets')));
            
            % load an rEFH
            params = setParams('datatype','spikecounts',...
                'mods',mods(iSession),'datafile',[localdir,filename],...
                'Nmsperbin',this_binwidth,...
                'trainingtime',this_trainingtime,...
                'fraction',this_fracneurons,...
                'hidsensoryratio',this_hidsensoryratio,...
                'phidtarget',this_phidtarget,...
                'SensoryUnitType',sensory_unit_type);
            
            % load the LTI system acquired with EM--or ask for a new one
            if ~TRAINNEWLDS&&any(ismember(decodernames,{'kfemstatic','kfemdynamic'}))
                Mstates = floor(params.numsUnits{1}/3);
                load(sprintf(...
                    '%sRBMish/BMI/EMparams/LDSOrd%03d_%s_%s_%imsBins_%03dsTrainingtime_%s',...
                    getdir('data'),Mstates,params.datatype,...
                    params.typeUnits{1}{1},...
                    params.Nmsperbin,...
                    params.trainingtime,filesuffix),...
                    'LDSparamsEM');
            else
                LDSparamsEM = [];   Bxz_EM = [];
            end
            
            % load the trained rEFH--or ask for a new one
            if ~TRAINNEWREFHS&&any(ismember(decodernames,{'refhstatic',...
                    'refhdynamic','refhstatic_stdnrml','refhdynamic_stdnrml'}))
                load(sprintf(...
                    '%sRBMish/EFHs/spikecounts/wts_%s_%s_%iHid_%imsBins_%03dsTrainingtime_%s',...
                    getdir('data'),params.datatype,...
                    params.typeUnits{1}{1},...
                    params.numsUnits{2},...
                    params.Nmsperbin,...
                    params.trainingtime,...
                    filesuffix),'wts','params');
                
                if strcmp(dataclass,'gpuArray')
                    wts{1} = gpuArray(wts{1}); wts{2} = gpuArray(wts{2});
                end
            else
                wts = [];
            end
            
            % run the filters
            [Rsqs(iSession,:,iSweptparam,:), NdataTest(iSession,1,iSweptparam)] =...
                filtersForNeuralData(wts,params,LDSparamsEM,decodernames,TOPLOT);
         
            % save the results
            out_file = sprintf('Rsqs_%s_%s_%s_%s.mat',...
                sweep_name,sensory_unit_type,monkey,date);
            Rsqs = gather(Rsqs); NdataTest = gather(NdataTest);
            save(out_file, 'Rsqs','decodernames','dates','seqnums',...
                'NdataTest','binwidths','trainingtimes','fracneurons',...
                'hidsensoryratios','phidtargets','monkey');
            %%%Nneurons(iSession,1,iBinsize) = params.numsUnits{1}(1);
            %%% why not also save kinemat and Nneurons??
            if strcmp(dataclass, 'gpuArray')
                Rsqs = gpuArray(Rsqs); NdataTest = gpuArray(NdataTest);
            end
            
        end
        try
            sendmail('makin@phy.ucsf.edu',...
                sprintf('finished training session %i of %s',...
                iSession,out_file));
        catch ME
            fprintf('failed to send e-mail...\n');
        end
    end
end








