% DBN   deep belief net training
%   DBN creates (see setParams.m) and trains a deep belief net.
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
%   -from Hinton and Salakhutdinov
%   by JGM
%-------------------------------------------------------------------------%

% init
clear all; close all

% params
params = setParams;
[~,machine] = system('hostname');
params.machine = strtrim(machine);
numsUnits = params.numsUnits;
numRBMs = length(numsUnits)-1;
Ncases = params.Ncases;
Nbatches = 40000/params.Ncases;
Nvis = numsUnits(1);
datagenargs = {};
TESTDECODING = 1;
Ntest = 5;

% say
paramDisplay(params);

%%% you need this line if you're doing backprop
% [testbatchdata0,xtest] = DATAGENPP(500,params,'correlation',0.8);

% pre-allocate wts cell array
switch params.nettype
    case 'ENCODER', wts = cell(numRBMs*2,1);
    case 'CLASSIFIER', wts = cell(numRBMs,1);
    otherwise, error('unrecognized type of network (see params.nettype)');
end
    
% pretraining
for i_rbm = 1:numRBMs
    Nhid = numsUnits(i_rbm+1);
    fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',i_rbm,Nvis,Nhid);
    restart = 1;
    
    % train
    tic; EFH; toc;
    
    % pack together weights for saving (hid => recog., vis => gener.)
    switch params.nettype
        case 'ENCODER'
            wts{i_rbm} = [vishid; hidbiases];
            wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases'];
            filename = 'EncoderWtsFile';
            save(filename,'i_rbm','numsUnits','wts','params','epoch');
        case 'CLASSIFIER'
            wts{i_rbm} = [vishid; hidbiases];
            save ClassifWtsFile i_rbm numsUnits wts params
    end

    % for next time through
    Nvis = Nhid;
    batchdata = batchphidmeans;   
    
end

% fine-tuning
% switch params.nettype
%     case 'ENCODER'
%         backprop;
%     case 'CLASSIFIER'
%         backpropclassify
% end
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
