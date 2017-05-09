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
%decodernames = {'static','kfobs','kfemstatic','kfemdynamic','ukf',...
%    'refhstatic','refhdynamic'};
decodernames = {'refhstatic','refhdynamic'};



TOPLOT = 0;
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
monkey = 'Indy';
Nmsperbin = 64; %%%[16,32,64,128];
Nbinsizes = length(Nmsperbin);
Ndecoders = length(decodernames);
Nstates = 6;
Nmspers = 1000; % fact
TRAINNEWREFHS = 1;
TRAINNEWLDS = 1;




% params
Nreachesmin = 300;
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
Rsqs = zeros(Nsessions,Nstates,Nbinsizes,Ndecoders,dataclass);
NdataTest = zeros(Nsessions,1,Nbinsizes,dataclass);


% loop across sessions and bins
for iSession = 1:Nsessions
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
    serverloc = sprintf('../../sabes_data2/joeyo/data/%s/%s/%s/%s/processed/',...
        lower(monkey),year,month,session);
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
    
    
    for iBinsize = 1:Nbinsizes
        
        % load an rEFH
        params = setParams('datatype','spikecounts',...
            'mods',mods(iSession),'datafile',[localdir,filename],...
            'Nmsperbin',Nmsperbin(iBinsize),'trainingtime',320,....
            'fraction',1,...
            'SensoryUnitType','Poisson');
            
      
        % load the LTI system acquired with EM--or ask for a new one
        if TRAINNEWLDS
            LDSparamsEM = [];   Bxz_EM = [];
        else
            Mstates = floor(params.numsUnits{1}/3);
                load(sprintf('%sRBMish/BMI/LDSOrd%03d_%s_%s_%imsBins_%03dsTrainingtime_%s',...
                getdir('data'),Mstates,params.datatype,params.typeUnits{1}{1},...
                Nmsperbin(iBinsize),params.trainingtime,filesuffix),...
                'LDSparamsEM');
        end
        
        % load the trained rEFH--or ask for a new one
        if TRAINNEWREFHS
            wts = [];
        else
            load(sprintf('%sRBMish/EFHs/spikecounts/wts_%s_%s_%iHid_%imsBins_%03dsTrainingtime_%s',...
                getdir('data'),params.datatype,params.typeUnits{1}{1},...
                params.numsUnits{2},Nmsperbin(iBinsize),...
                params.trainingtime,filesuffix),'wts','params');
            if strcmp(dataclass,'gpuArray')
                wts{1} = gpuArray(wts{1}); wts{2} = gpuArray(wts{2});
            end
        end
        
%         % run the filters
%         [Rsqs(iSession,:,iBinsize,:), NdataTest(iSession,1,iBinsize)] =...
%             filtersForNeuralData(wts,params,LDSparamsEM,decodernames,TOPLOT);
%         %%%Nneurons(iSession,1,iBinsize) = params.numsUnits{1}(1);
%         
%         % save results
%         save(sprintf('allRsqs%s%03dsTrainingtime_%s_%s',...
%             monkey,params.trainingtime,date,params.typeUnits{1}{1}),...
%             'decodernames','Rsqs','dates','seqnums','NdataTest');
%         %%% kinemat?  trainingTime? binwidths? etc.
        params.mods
    end

    
    %%% forget about deleting them afterwards....
    
end


