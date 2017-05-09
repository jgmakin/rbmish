% nonflatpriorCond.m
%   Compute the *conditional* error statistics for a few different
%   directions (see conditionals.m) under the Gaussian-prior-trained model.
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%   -re-created nonflatCondStats.mat from wtsPrior050.mat
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


clear; clc; close all

% load
% wtsfile = 'results/wtsPrior.mat';
wtsfile = 'results/finalwts/wtsPrior050.mat';
% NNfile = 'results/NNPrior.mat';
condfile = 'results/finalwts/nonflatCondStats';

% test
RESTART = 0;
inc = 0.01;
conditionals