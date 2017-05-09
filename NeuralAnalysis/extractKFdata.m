function extractKFdata(PP_DATA,NEURAL,All_TargPts,tag)
% load and collect the data you'll need
%
% USAGE: 
%{   
   tag = 'D080616'

   datapath = 'C:\#DATA\dmitri_datafiles\';
   fileprefix = [datapath,tag];
 
   load([fileprefix,'_neural.mat']);
   load([fileprefix,'_bparams.mat']);
   load([fileprefix,'_behav.mat']);
 
   extractKFdata(PP_DATA,NEURAL,All_TargPts,tag);
%}
   
%-------------------------------------------------------------------------%
% Revised: 06/04/13
%   -gave it some arguments.  Now you call it a la scrap.m (see usage note)
% Cribbed: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%

% tuning trials
TunTrial = cat(1,PP_DATA(:).trialtype)==1;
Nttrials = sum(TunTrial);

[St(1:Nttrials).pos] = deal(PP_DATA(TunTrial).Reach_svect_traj);
[St(1:Nttrials).vel] = deal(PP_DATA(TunTrial).Reach_svect_vel);
[St(1:Nttrials).acc] = deal(PP_DATA(TunTrial).Reach_svect_acc);
[St(1:Nttrials).t] = deal(PP_DATA(TunTrial).Reach_traj_time);

[UnitSpikesT(1:length(NEURAL)).t] = NEURAL(:).UnitSpikes;
[UnitSpikesT(1:length(NEURAL)).id] = NEURAL(:).unitID;

targetLocsT = unique(All_TargPts,'rows');

filename = ['KFtuningdataHHS',tag];
save(filename,'St','UnitSpikesT','targetLocsT')


% reaching trials
ReachTrial = cat(1,PP_DATA(:).trialtype)==2;
Nrtrials = sum(ReachTrial);

[Sr(1:Nrtrials).pos] = deal(PP_DATA(ReachTrial).Reach_svect_traj);
[Sr(1:Nrtrials).vel] = deal(PP_DATA(ReachTrial).Reach_svect_vel);
[Sr(1:Nrtrials).acc] = deal(PP_DATA(ReachTrial).Reach_svect_acc);
[Sr(1:Nrtrials).t] = deal(PP_DATA(ReachTrial).Reach_traj_time);

[UnitSpikesR(1:length(NEURAL)).t] = NEURAL(:).UnitSpikes;
[UnitSpikesR(1:length(NEURAL)).id] = NEURAL(:).unitID;

targetLocsR = unique(All_TargPts,'rows');

filename = ['KFreachingdataHHS',tag];
save(filename,'Sr','UnitSpikesR','targetLocsR')


end

