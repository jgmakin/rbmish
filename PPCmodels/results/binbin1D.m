% binbin1D
%-------------------------------------------------------------------------%
% (3) 1D binomial-binomial model with nonlinearity (cosine)
%
% This is not perfect---which is strange, b/c the 2D version is pretty darn
% good.  There's still some nonzero bias.
%-------------------------------------------------------------------------%

clear; clc; close all

% load
load results/new/wtsBB1D % results/wts1DBB
% load results/NN1DBB

% test
% [ErrorStats net] = nnDecode(wts,'prop',params,'pretrained',net);
ErrorStats = test(wts,params);