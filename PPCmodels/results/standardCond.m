clear; clc; close all

% load
wtsfile = 'results/new/wtsPrior.mat';
condfile = 'results/new/sixCondStatsPrior8std';
% condfile = 'results/new/sixCondStatsPrior4std';

% test
RESTART = 1;
inc = 8/150; % 0.5; % this will be overwritten if RESTART = 0;
conditionals