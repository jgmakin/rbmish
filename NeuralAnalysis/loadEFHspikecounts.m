function [model,params,model_full_path] = loadEFHspikecounts(monkey,iSession,...
    binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,model_type)
% loadEFHspikecounts    Load wts and params for a saved 'spikecounts' model
%
% USAGE:
%   [wts,params] = loadEFHspikecounts(monkey,iSession,...
%       binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb);
%
% You have run and saved so many 'spikecounts' models that it's become 
% useful to have a way of converting the parameters into the corresponding
% filename (and then loading the file).

%-------------------------------------------------------------------------%
% Created: 10/20/17
%   by JGM
%-------------------------------------------------------------------------%

% load data file
results_filename = sprintf('%sRBMish%sBMI%sRsqs_%sSweep_%s_%s.mat',...
    getdir('data'),filesep,filesep,...
    [upper(swept_param(1)),swept_param(2:end-1)],...
    obsv_dstrb,monkey);
load(results_filename,'dates')
session = dates{iSession};
switch monkey
    case 'Chewie'
        year    = session(1:4);
        month   = session(5:6);
        day     = session(7:8);
        datafile = sprintf('%s_%s%s%s.mat',monkey,month,day,year);
    otherwise
        load(results_filename,'seqnums') % you need this, too
        date_seqnum = sprintf('%s_%02i',session,seqnums(iSession));
        datafile = sprintf('spikes_and_kinematics_%s.mat',date_seqnum);
end
datadir = sprintf('%s%s_datafiles%s',getdir('data'),monkey,filesep);
if ~exist([datadir, datafile],'file')
    fprintf('data file not found, returning\n');
    return
end

switch model_type
    case 'EFH'
        % for all decoders, you need the params from the refh file
        model_dir = sprintf('%sRBMish%sEFHs%sspikecounts%s',...
            getdir('data'),filesep,filesep,filesep);
        model_filename = sprintf('wts_spikecounts_%s_*Hid_%imsBins_%isTrainingtime_%s.mat',...
            obsv_dstrb, binwidth, trainingtime, date_seqnum);
        model_filename = filenameForFracneurons(model_dir,model_filename,...
            fracneuron,obsv_dstrb,monkey);
        model_full_path = [model_dir, model_filename];
        
    case 'LDS'
        model_dir = sprintf('%sRBMish%sBMI%sEMparams%s',...
            getdir('data'),filesep,filesep,filesep);
        model_filename = sprintf('LDSOrd*_spikecounts_%s_%imsBins_%isTrainingTime_%s.mat',...
            obsv_dstrb,binwidth,trainingtime,date_seqnum);
        model_filename = filenameForFracneurons(model_dir,model_filename,...
            fracneuron,obsv_dstrb,monkey);
        model_full_path = [model_dir, model_filename];
        
    otherwise
        error('unrecognized BMI model -- jgm');
end
        
% now (try to) load the parameters from this file
if exist(model_full_path,'file')~=2
    fprintf('data file not found, returning...\n');
    model = [];
else
    load(model_full_path);
    switch model_type
        case 'EFH'
            model = wts;
        case 'LDS'
            model = LDSparamsEM;
    end
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function model_filename = filenameForFracneurons(model_dir,...
    model_filename,fracneuron,obsv_dstrb,monkey)

% the training times are reflected in the number of hidden units/states
matching_filenames = ls([model_dir, model_filename]);
Nfiles = size(matching_filenames,1);

% find the number of hiddens in all the otherwise matching file names
numhid = zeros(Nfiles,1);
for iFile = 1:Nfiles
    this_file = matching_filenames(iFile,:);
    switch this_file(1:3)
        case 'wts'
            i0 = length(sprintf('wts_spikecounts_%s_',obsv_dstrb)) + 1;
            iF = strfind(this_file,'Hid') - 1;
        case 'LDS'
            i0 = 7;
            iF = 9;
        otherwise
            fprintf('WARNING: unexpectedly found a file name ')
            fprintf(' without LDS or wts\n');
    end
    numhid(iFile) = str2double(this_file(i0:iF));
end
switch Nfiles
    case 4 % all trainingtimes files are present
        
        % load the canonical fracneurons sweep
        load(sprintf('%sRBMish%sBMI%sRsqs_FracneuronSweep_%s_%s.mat',...
            getdir('data'),filesep,filesep,obsv_dstrb,monkey),'fracneurons')
        
        % sort Nhiddens, grab file that matches fracneuron
        [~,sorted_inds] = sort(numhid);
        iFile = sorted_inds(fracneurons==fracneuron);  
        model_filename = matching_filenames(iFile,:);
        
    case 1 % probably downloaded for this purpose
        fprintf('WARNING: only found refh wts for one fracneuron.')
        fprintf('  Assuming this the correct file and proceeding...\n')
        model_filename = matching_filenames(Nfiles,:);
        
    otherwise
        fprintf('%i possible matches (out of four); don''t know',Nfiles)
        fprintf(' which fracneuron is intended so exiting\n')
        model_filename = [];
        return
end


end
%-------------------------------------------------------------------------%
