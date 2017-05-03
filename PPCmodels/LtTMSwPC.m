% LtTMSwPC
%
% Reproduce figures from Makin2015b:
% "Learning to Estimate Dynamical State w/Probabilistic Population Codes"
%
% This is an ongoing project and not finished.  Ultimately it will replace
% (or call) many of the m-files in the RBM/results folder.  The goal is to
% generate tikz figures, integrated with your color and symbol tex schemes,
% especially rbmish.sty.


%-------------------------------------------------------------------------%
% Revised: 03/04/16
% Created: 02/29/16
%   -by JGM
%-------------------------------------------------------------------------%


%% retrain with recurrentEFHtrain.m.......






%% Figure 1D:
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_ManyXprmts.mat'
plotIllustrativeTrajectories('LTI-PPC','rEFH',1,{})
% figure saved as
%   illustrativeTraj_1DrEFH_Joint-Angle_[DATE].tex



%% Figure 1E:
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_ManyXprmts.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_ManyXprmts.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_ManyXprmts.mat'
filterNames = {'sensory','LDSOrd1','rEFH','LDSOrd2','KFtrue'};
LDStest('rEFH',1,{},'Xprmts',filterNames);
% figure saved as
%   'MSEsBar_Joint-Angle_[DATE].tex





%% Figure 2: Mean squared errors for various dynamical systems
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_ManySprings.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_ManySprings.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_ManySprings.mat'
filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
LDStest('rEFH',1,{},'Springs',filterNames);
% Figure 2A is saved as 
%   'MSEsVsSprings_Joint-Angle_[date].tex'

% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_ManyDampers.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_ManyDampers.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_ManyDampers.mat'
filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
LDStest('rEFH',1,{},'Dampers',filterNames);
% Figure 2B is saved as 
%   'MSEsVsDampers_Joint-Angle_[date].tex'

% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_ManyMasses.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_ManyMasses.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_ManyMasses.mat'
filterNames = {'LDSOrd1','rEFH','LDSOrd2','KFtrue'};
LDStest('rEFH',1,{},'Masses',filterNames);
% Figure 2C is saved as 
%   MSEsVsMasses_Joint-Angle_[date].tex





%% Figures 3C,E
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_1DrEFHwithEC_ManyXprmts.mat'
plotIllustrativeTrajectories('LTI-PPC','rEFH',1,{'withEC'})





%% Figures 3D,F
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_1DrEFHwithEC_ManyXprmts.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd1_1DrEFHwithEC_ManyXprmts.mat'
%   '[dataDir/]RBMish/EMparams/LDSOrd2_1DrEFHwithEC_ManyXprmts.mat'
filterNames = {'sensory','KFobsNoCtrl','LDSOrd2','rEFH','LDSOrd3','KFtrue'};
LDStest('rEFH',1,{'withEC'},'Xprmts',filterNames);
% Figure 3D saved as:
%   'MSEsBar_Joint-Angle_[DATE].tex
% Figure 3F saved as:
%   'MSEsBar_Efference-Copy_[DATE].tex




%% Figure 4: Box-and-whisker plot of MSEs for EM-based models and rEFHs
% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_VaryingHids.mat'
%	'[dataDir/]RBMish/EFHs/wts_TRBM_1D_LTI-PPC_VaryingHids.mat'
%	'[dataDir/]RBMish/EFHs/wts_RTRBM_1D_LTI-PPC_VaryingHids.mat'
% 	'[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_VaryingHids.mat'
% 	'[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_VaryingHids.mat'
%	'[dataDir/]RBMish/testdata/LDSdata_1D_LTI-PPC.mat'
boxplotREFH('LTI-PPC',1,{});

% Requires the following saved files:
%	'[dataDir/]RBMish/EFHs/wts_rEFH_1D_LTI-PPC_withEC_VaryingHids.mat'
% 	'[dataDir/]RBMish/EMparams/LDSOrd1_1D_LTI-PPC_withEC_VaryingHids.mat'
% 	'[dataDir/]RBMish/EMparams/LDSOrd2_1D_LTI-PPC_withEC_VaryingHids.mat'
%	'[dataDir/]RBMish/testdata/LDSdata_1D_LTI-PPC_withEC.mat'
boxplotREFH('LTI-PPC',1,{'withEC'});

% Figure 4a is saved as 1D-LTI-PPC_boxplotsB_[date].tex 
% Figure 4b is saved as 1D-LTI-PPC_withEC_boxplotsB_[date].tex 
% REFH_JMLR Figure 4 is saved as 1DrEFH_ErrorStatsVsNumHiddens_[date].tex 




%% Figure 5: Position and velocity RFs of hidden units.

%%%% this will involve wts_1DrEFH_NoSpring_ManySpeeds.mat































