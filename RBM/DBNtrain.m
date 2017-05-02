% DBNtrain   deep belief net training,'
%   DBNtrain creates (see setParams.m) and trains a deep belief net.

%-------------------------------------------------------------------------%
% Revised: 09/30/10
%   -cleaned up
% Revised: 05/25/10
%   -made number of layers variable
%   -consolidated rbmhidlinear with rbm (automatically makes the "top"
%       layer linear)
% Revised: 05/24/10
%   -fixed formating
% Adapted: 05/24/10
%   -from Hinton and Salakhutdinov (see bottom of file)
%   by JGM
%-------------------------------------------------------------------------%

% init
clear; % close all

% params
% params = setParams('datatype','2Dinteg');
% params = setParams('datatype','bouncingballs','algorithm','RTRBM');
params = setParams('datatype','spikecounts','mods',{'M1S1'},...
    'datafile','Jackson_datafiles/spikes_and_kinematics_jackson_A_45.mat',...
    'Nmsperbin',64,'trainingtime',320,'SensoryUnitType','Poisson');

if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

numsUnits = params.numsUnits;
datagenargs = {};
if isfield(params,'EACHBATCHISATRAJ')
    EACHBATCHISATRAJ = params.EACHBATCHISATRAJ;
else
    EACHBATCHISATRAJ = 0; 
end
if isfield(params,'Npretrain')
    Npretrain = params.Npretrain;
else
    Npretrain = 0;
end
if isfield(params,'sparse')
    sparsitycost= params.sparse.cost;    
    estimRate   = params.sparse.phidNewFrac;
    phidTarget  = params.sparse.phidTarget;
else
    sparsitycost= 0;
end
getLatents  = params.getLatents;
getData     = params.getData;
NEFHs       = length(numsUnits)-1;
Ncases      = params.Ncases;
Nbatches    = params.Nbatches;
Ntest       = params.NepochsMax + 1;
%%%Ntest       = 5;
NepochsMax  = params.NepochsMax;



% init
paramDisplay(params);
wts = cell(NEFHs*2,1);

    
% pretraining
for iEFH = 1:NEFHs
    fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',...
        iEFH,sum(numsUnits{iEFH}),sum(numsUnits{iEFH+1}));
    RESTART = 1;
    
    % train
    tic; EFH; toc;
    
    % pack together weights for saving (hid => recog., vis => gener.)
    wts{iEFH} = [vishid; hidbiases];
    wts{NEFHs*2-iEFH+1} = [vishid'; visbiases'];
    filename = 'EncoderWtsFile';
    save(filename,'iEFH','numsUnits','wts','params','epoch');
    
end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton  
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with a
% note saying that the original programs are available from our web page. 
% The programs and documents are distributed without any warranty, express
% or implied.  As the programs were written for research purposes only,
% they have not been tested to the degree that would be advisable in any
% important application.  All use of these programs is entirely at the
% user's own risk.
%
%
% This program pretrains a deep autoencoder for MNIST dataset You can set
% the maximum number of epochs for pretraining each layer and you can set
% the architecture of the multilayer net.
%-------------------------------------------------------------------------%

