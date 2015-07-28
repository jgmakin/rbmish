% DBNRESTART    Adds a deeper layer to an existing DBN and trains it
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

% path
path(path,'../ee125/fxns');
path(path,'parallel');
path(path,'results');
path(path,'dynamical');
path(path,'tuningcurves');
path(path,'retired');
path(path,'scratch');
path('../utils',path);                      % contains tex.m
path(path,'../tools')

% update dbn params
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 you need to set these
newlayerlength = params.numsUnits(1);  % 3200; % size(wts{1},2);
newlayertype = 'BP';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.swing = 0.2;                          % just in case
maxepoch = params.DBNmaxepoch;
params.numsUnits = [params.numsUnits newlayerlength];
params.typeUnits = [params.typeUnits newlayertype];
numsUnits = params.numsUnits;
numRBMs = length(numsUnits)-1;
i_rbm = numRBMs;
numvis = params.numsUnits(end-1);

% enlarge wts (cell) array
W = wts;
wts = cell(numRBMs*2,1);
wts([1:numRBMs-1,(numRBMs+2):end]) = W([1:numRBMs-1,numRBMs:end]);
clear W;

% make data
datagenargs = {}; %%%% but you may want to add stuff here!
datagenargs = [datagenargs,{'dbndepth',i_rbm,'dbnwts',wts}];
[batchdata,Sbatch] = DATAGENPP(1000,params,datagenargs{:});
[Ncases,Nvis,Nbatches] = size(batchdata);
Nvis = params.numsUnits(end-1);
Nhid = params.numsUnits(end);

% run
restart = 1;
EFH


















