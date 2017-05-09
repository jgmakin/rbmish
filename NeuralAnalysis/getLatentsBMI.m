function [X,Q,badneurons] = getLatentsBMI(Nexamples,dataclass,T,testortrain,params,varargin)
% getLatentsBMI
%
% USAGES:
%
%   [X,Q] = getLatentsBMI(Nexamples,yrclass,T,'test',params);
%
%   [X,Q] = getLatentsBMI(Nexamples,yrclass,T,'train',params,...
%       'sequencelength',1000);
%
%   [X,Q,badneurons] = getlatentsBMI(...)
%

%-------------------------------------------------------------------------%
% Revised: 05/08/17
%   -returns badneurons, the neurons which are culled
% Revised: 04/17/17
%   -added code to remove units that happen to have identical firing
%   patterns ("duplicates") over either the training or testing portions.
% Revised: 04/13/17
%   -incorporated test-data trimming into this file
%   -changed the way single sequences are requested
%   -moved the normalization for StandardNormal units into here (see note
%   below)
% Created: 04/21/16
%   by JGM
%-------------------------------------------------------------------------%

%%%% TO DO:
% (1) At the moment, this doesn't work for HHS data.  See the m-files in
% the directory NeuralAnalysis for how to deal with those....

% necessary parameters
T = defaulter('sequencelength',T,varargin{:});
trainingtime = params.trainingtime;
mods = params.mods;
datafile = params.datafile;
Fs = params.Fs;
Nmsperbin = params.Nmsperbin;
Nmspers = 1000; % fact
sperbin = (Nmsperbin/Nmspers);
Nsamplesperbin = round(Fs*sperbin);

% load neural data
[spikedata, St, Ndims] = loadNeuralData(dataclass,datafile,mods,Fs);
minRate = 0.5;                              % Hz.  Somewhat arbitrary

% useful params
BinParams.m = Nsamplesperbin;
BinParams.dt = 1/Fs;
BinParams.Ndims = Ndims;
BinParams.Nstates = size(St(1).X,2);
BinParams.BINMETHOD = 'void';

% assemble kinematic state
if 0
    kinematicsFilterOrder = 2;
    kinematicsFilterCutoff = 6;
    [bbb,aaa] = butter(kinematicsFilterOrder,kinematicsFilterCutoff/(Fs/2));
    X = filtfilt(bbb,aaa,X);
end

% assemble usable "neurons"
UnitSpikesT = cell2struct(spikedata(:), 't', 2);

% transform kinematic data and spike times into useful data
[R,X,~] = binSpikeCounts(St,UnitSpikesT,BinParams);

too_slow = mean(R)/sperbin <= minRate;
dead_neurons = arrayfun(@(ii)(isempty(spikedata{ii})),1:numel(spikedata));
badneurons = too_slow | dead_neurons;

R(:,badneurons) = [];

if isfield(params,'fraction'), R = R(:,1:(floor(end*params.fraction))); end


% find the split point for training and testing data
switch trainingtime
    case 'half'
        NbinsTrain = floor(size(R,1)/2);
    otherwise % training time is a number
        NbinsTrain = floor(trainingtime/sperbin);
end



% which bins? which neurons?
binsTrain   = 1:NbinsTrain;
binsTest    = getTestBins(NbinsTrain,size(X,1),sperbin,params.datafile);

% for numerical stability, find neurons that happen to duplicate
trainDuplicates = findDuplicates(R(binsTrain,:));
testDuplicates  = findDuplicates(R(binsTest,:));
duplicates = [trainDuplicates; testDuplicates];
fprintf('removing %i ''duplicate'' neuron(s)...\n',length(duplicates));
neurons = logical(prod((1:size(R,2))~=duplicates,1));


% now sub-select bins, neurons
switch testortrain
    case 'train', bins = binsTrain;
    case 'test',  bins = binsTest;
    otherwise
        fprintf('invalid data request! -- returning *all* data\n');
        bins = 1:size(X,1);
end
R = R(bins,neurons);
X = X(bins,:);


% params.getData may not have been created yet
switch params.typeUnits{1}{1} %%% assume only one
    case 'StandardNormal'
        % Where should this happen?  It would be "cheating" to normalize
        % the training and testing data together, since that allows testing
        % info to leak into the training.  On the other hand, normalization
        % should be applied to the whole trajectory, not the random set of
        % minibatches assembled below.
        R = (R - mean(R))./std(R);
        %%% R = (R - mean(R))/chol(cov(R));
    case 'Bernoulli'
        R(R>1) = 1;
end


% now select out Ntraj trajectories of length T
if strcmp(T,'singlesequence'), T = size(X,1); Nexamples = size(X,1); end
Ntraj       = floor(Nexamples/T);
NbinsTotal  = size(X,1);
i0          = ceil((NbinsTotal-T)*rand(Ntraj,1));
inds        = (i0 + (1:T))';
X           = X(inds,:);
Q.restarts  = 1:T:(Ntraj*T);
Q.R         = R(inds,:);
Q.T         = T;


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [spikes,St,Ndims] = loadNeuralData(yrclass,datafile,mods,Fs)


if any(strcmp(datafile(1:4),{'Indy','Loco','Jack'}))
    % chan_names        cell array of names, ie 'M1 001'
    % cursor_pos        cursor position, sampled at 250 Hz [t by 2]
    % finger_pos        finger position, 250 Hz [t by 3]
    % spikes            cell array of spike times (seconds) {channel,unit}
    % t                 timebase for kinematics, seconds

    % hard-coded
    Norder = 3;
    
    % load
    load([getdir('data'),datafile],'cursor_pos','spikes','t');
    X = cursor_pos;
    Ndims = size(X,2);
    
    % just one "trial"---may be different for other data, e.g. HHS
    St(1).X = X;
    for iOrder = 2:Norder
        X = diff(X)*Fs;
        X = [X; X(end,:)];
        St(1).X = [St(1).X,X];
    end
    St(1).t = t;
    
    
elseif strcmp(datafile(1:6),'Chewie')
    
    % load
    load([getdir('data'),datafile],'out_struct');
    
    % spikes
    IDs = cat(1,out_struct.units(:).id);
    Nelectrodes = max(IDs(:,1));
    Nunits = length(out_struct.units);
    spikes = cell(Nelectrodes,1); %%% but there are probably more cols
    for iUnit = 1:Nunits
        thisUnit = out_struct.units(iUnit);
        jElectrode = thisUnit.id(1);
        jUnit = thisUnit.id(2) + 1; %%% start w/0
        if jUnit ~= 256
            spikes{jElectrode,jUnit} = thisUnit.ts;
        end        
    end
    
    % kinematics
    St(1).X = cat(2,out_struct.pos(:,2:3),out_struct.vel(:,2:3),...
        out_struct.acc(:,2:3));
    St(1).t = out_struct.pos(:,1);
    
    % number of spatial dimensions
    Ndims = 2;
    
else
    error('unexpected BMI data file -- jgm');
end

% GPU?
if strcmp(yrclass,'gpuArray')
    St(1).X = gpuArray(St(1).X);
    for i=1:size(spikes,1)
        for j=1:size(spikes,2)
            spikes{i,j} = gpuArray(spikes{i,j});
        end
    end
    St(1).t = gpuArray(St(1).t);
end



%%% might want to re-write this so that using both M1 and S1 data is
%%% labeled as two "modalities," rather than one with the name 'M1S1'.
details = whos(matfile([getdir('data'),datafile]),'chan_names');
switch mods{1} %%% assume only one modality....
    case {'M1','S1'}
        if isempty(details)
            fprintf('warning: no channel names; assuming all M1\n');
        else
            load([getdir('data'),datafile],'chan_names');
            spikes = spikes(strncmp(chan_names,mods,2),:);
        end
        
    otherwise
        % do nothing
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function bins = getTestBins(NbinsTrain,binF,sperbin,datafile)
% Take in *all* the data, and return the (trimmed) test data

% ceteris paribus, test data start right after the training data end
bin0 = NbinsTrain + 1;

% get the region to trim (these were found "by hand")
switch datafile
    case 'Chewie_datafiles/Chewie_10032013.mat'
        fprintf('trimming unfortunate Kording data...\n');
        badsecond0 = 508.8;
        bin0 = ceil(badsecond0/sperbin);
    case 'Chewie_datafiles/Chewie_12192013.mat'
        fprintf('trimming unfortunate Kording data...\n');
        badsecondF = 620.8;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170210_03.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 480;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170213_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 1996.8;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170214_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 695;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170215_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecond0 = 346.24;
        bin0 = ceil(badsecond0/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170216_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 770;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170217_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecond0 = 337.4;
        badsecondF = 720;
        bin0 = ceil(badsecond0/sperbin);
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170227_04.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 1965;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170228_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 1285;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170301_05.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 590;
        binF = floor(badsecondF/sperbin);
    case 'Loco_datafiles/spikes_and_kinematics_20170302_02.mat'
        fprintf('trimming unfortunate Loco data...\n');
        badsecondF = 1450;
        binF = floor(badsecondF/sperbin);
end

% finally, trim
bins = bin0:binF;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function iDuplicates = findDuplicates(R)

try
    [ix,iy] = find(squeeze(sum(R == permute(R,[1,3,2]),1)) == size(R,1));
    iDuplicates = iy(ix<iy);   % just from the lower triangle
catch ME
    fprintf('not enough memory to check for ''duplicate'' neurons--skipping\n');
    iDuplicates = zeros(0,1);
end

end
%-------------------------------------------------------------------------%
