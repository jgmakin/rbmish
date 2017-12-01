% DBNrestart    Adds a deeper layer to an existing DBN and trains it
%-------------------------------------------------------------------------%
% Revised: 05/14/13
%   -changed to accommodate the new DATAMAKER---but didn't check
% Revised: 09/30/10
%   -changed training data to variable gain
% Created: 09/15/10
%   by JGM
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%
% load wtsPartial
% or 
% wts = wts1;
% params = params1;
%%%%%%%%%%%%%


% update dbn params
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 you need to set these
newlayerlength = params.numsUnits{2};
newlayertype = {'Bernoulli'};
params.Ncdsteps = 1;
params.NepochsMax = 90;
params = setLearningSchedules(50,12,12,'exp',params);
datagenargs = {}; %%%% but you may want to add stuff here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
NepochsMax = params.NepochsMax;
params.numsUnits = [params.numsUnits {newlayerlength}];
params.typeUnits = [params.typeUnits {newlayertype}];
numsUnits   = params.numsUnits;
numRBMs     = length(numsUnits)-1;
iRBM        = numRBMs;
EACHBATCHISATRAJ = 0;
Ntest       = 5;
Npretrain   = 0;
Ncases      = params.Ncases;
Nbatches    = params.Nbatches;


% enlarge wts (cell) array
W = wts;
wts = cell(numRBMs*2,1);
wts([1:numRBMs-1,(numRBMs+2):end]) = W([1:numRBMs-1,numRBMs:end]);
clear W;

% make data
datagenargs = [datagenargs,{'dbndepth',iRBM,'dbnwts',wts}];
batchdata = generateData(Ncases*Nbatches,params,datagenargs{:});

% run
restart = 1;
EFH

% store weights
wts{iRBM} = [vishid; hidbiases];
wts{numRBMs*2-iRBM+1} = [vishid'; visbiases'];
filename = 'EncoderWtsFile';
save(filename,'iRBM','numsUnits','wts','params','epoch');















