clear gv1;
iTraj = ceil(rand*size(LDSdataTest.Z,1)); 
% run number; arbitrary number from 1 to 40

states = squeeze(LDSdataTest.Z(iTraj,:,:));
input_act = squeeze(LDSdataTest.R(iTraj,:,:));
%%%%
V0 = squeeze(LDSdataTest.R);
%%%%
hidden = squeeze(V0(iTraj,:,:));

gv1 = grid_viewer2(states,params,input_act,hidden);
gv1.play(10)