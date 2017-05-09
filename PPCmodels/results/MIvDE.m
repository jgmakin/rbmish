% MIvDE
%
% Reproduce figures from Makin2013b:
% "Learning Multisensory Integration and Coordinate Transformation via
% Density Estimation"
%
% This is an ongoing project and not finished.  Ultimately it will replace
% (or call) many of the m-files in the RBM/results folder.  The goal is to
% generate tikz figures, integrated with your color and symbol tex schemes,
% especially rbmish.sty.


%-------------------------------------------------------------------------%
% Revised: 02/29/16
% Created: 02/22/16
%   -by JGM
%-------------------------------------------------------------------------%


%% retrain models....
% multisensory integration
%   set params.datatype = '2Dinteg'
%       (1) run DBN
%       (2) non-flat prior: in DBN, add lines:
%           p0.cov = 1.0e-03*[0.2467,0;0,0.1097];
%           p0.mu = [-0.3927;1.5708];
%        and set datagenargs = {'stimulusprior',p0);};
%
% hierarchical integration
%   run hiertrain.m
%
% coordinate transformation
%   set params.datatype = '1Daddition', run DBN




%% Figure 1: schematic parts
% NB: not all the parts are here

clear; clc;
%%%load([getdir('data'),'RBMish/EFHs/wts_2Dinteg_160226.mat']);
load([getdir('data'),'RBMish/EFHs/wts_2Dinteg_900_31-Dec-2016.mat']);
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
r0 = generateData(1,params);
illustrateUnisensoryResponses(r0,params);
% figures saved as 
%   exPropResp-[date].tex
%   exVislResp-[date].tex




%% Figure 2: Recovery of the Posterior Mean
clear; clc;
load([getdir('data'),'RBMish/EFHs/wts_2Dinteg_160226.mat']);
tikzConditionalErrorEllipses(15,15,2.5,wts,params);
% figure saved as 
%   2DerrorStats-Joint-Angle-[date].tex

eStats = mastertest(wts,params);
% saved as 
%   2DconditionalErrorStats-Joint-Angle-[date].tex













%% plot the workspaces
clear; clc;
params = setParams('datatype','2Dinteg');
roboparams = params.roboparams;

L1 = roboparams.armlengths(1);
L2 = roboparams.armlengths(2);
thmin = roboparams.thmin;
thmax = roboparams.thmax;

M = 100;
inds = 1:M;
Th(inds,1) = linspace(thmin(1),thmax(1),M);
Th(inds,2) = thmin(2)*ones(M,1);

inds = inds + M;
Th(inds,1) = thmax(1)*ones(M,1);
Th(inds,2) = linspace(thmin(2),thmax(2),M);

inds = inds + M;
Th(inds,1) = linspace(thmax(1),thmin(1),M);
Th(inds,2) = thmax(2)*ones(M,1);

inds = inds + M;
Th(inds,1) = thmin(1)*ones(M,1);
Th(inds,2) = linspace(thmax(2),thmin(2),M);

X = FK2link(Th,roboparams,1);

figure(1); clf; hold on;
plot(X(:,1),X(:,2))
axis equal;

%%%
hold off;

figure(2); clf; hold on;
plot(Th(:,1),Th(:,2));
%axis equal;
hold off;

%%% might want to draw *the arm* in each of the "corners" of this plot!!!




































