function [pKF, LDSparams] = obsKFwithNoInput(LDSdataTrain,LDSdataTest,params)
% [pKF, LDSparams] = obsKFwithNoInput(LDSdata,LDSdataTest,params)
%
% Suppose the data are driven by an "efference copy," but the filter just
% decides to ignore it.  This function provides the parameters and filter
% cumulants in that case.

%-------------------------------------------------------------------------%
% Created: 09/04/14
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Ny,Nx] = size(params.dynamics.C);
Noutputs = size(LDSdataTrain.Y,2);

% train
LDSdataTrain.Z = LDSdataTrain.Z(:,1:Nx,:);
LDSdataTrain.Y = LDSdataTrain.Y(:,1:Ny,:);
LDSdataTrain.SigmaY = LDSdataTrain.SigmaY(:,1:Ny,1:Ny,:);
LDSparams = getLDSparams(params,'observed',LDSdataTrain);

% test
LDSdataTest.Z = LDSdataTest.Z(:,1:Nx,:);
LDSdataTest.Y = LDSdataTest.Y(:,1:Ny,:);
LDSdataTest.SigmaY = LDSdataTest.SigmaY(:,1:Ny,1:Ny,:);
pKF = KF4PPC(LDSdataTest,LDSparams,'obs');

keyboard

% pad with NaNs
[Ntraj,~,T] = size(pKF.Xpct);
pKF.Xpct = cat(2,pKF.Xpct,nan([Ntraj,Noutputs-Ny,T])); % 'like' ?...

end