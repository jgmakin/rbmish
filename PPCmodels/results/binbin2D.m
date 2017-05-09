% binbin2D
%-------------------------------------------------------------------------%
% (3) 2D binomial-binomial model
%   not as good as 1D, or as PB2D (standard), but still quite good
%-------------------------------------------------------------------------%

clear; clc; close all

% load
load results/wts2DBB
load results/NN2DBB

% test
[ErrorStats,net] = nnDecode(wts,params.NS,params,'pretrained',net);