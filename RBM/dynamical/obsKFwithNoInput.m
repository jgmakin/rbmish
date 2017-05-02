function [pKF, LDSparams] = obsKFwithNoInput(Ytrain,Xtrain,...
    unisensCmlntsTest,Q,params)
% USAGE:
%   [pKF, LDSparams] = obsKFwithNoInput(Ytrain,Xtrain,unisensCmlntsTest,...
%       Q,params)
%
% This function serves an extremely circumscribed purpose: fit and test an
% LTI system on only one of the observations and only the first two of the
% states.  You use it when the third (ignored) state is the "control," and
% the second observation is the "efference copy" of this control.  This is
% a pretty low bar, since there are in fact better ways to use just two
% states to fit these data--ways which EM will find, but this (regression-
% based) function will not.

%-------------------------------------------------------------------------%
% Revised: 01/09/17
%   -rewrote from scratch as part of Grand Revision
% Created: 09/04/14
%   by JGM
%-------------------------------------------------------------------------%

% Ns
T = Q.T;
Nexamples = size(Ytrain,1);
Ntraj = floor(Nexamples/T);
Ndims = size(unisensCmlntsTest.Xpct,2);

% use only the first modality
params.mods = params.mods(1);
unisensCmlntsTest.Xpct = unisensCmlntsTest.Xpct(:,:,1);
unisensCmlntsTest.Info = unisensCmlntsTest.Info(:,:,:,1);
Ytrain = Ytrain(:,1:size(unisensCmlntsTest.Xpct,2));

% use only the first two states
Xtrain = Xtrain(:,1:2);

% fit
LDSparams = learnfullyobservedLDS(shortdata(Ntraj,3,Ytrain),...
    shortdata(Ntraj,3,Xtrain));

% filter
pKF = KFposteriorization(unisensCmlntsTest,Q,LDSparams,params);

% pad with NaNs to make the same size as other posteriors
pKF.Xpct = cat(3,pKF.Xpct,NaN([Nexamples,Ndims,]));
pKF.Info = cat(4,pKF.Info,NaN([Nexamples,Ndims,Ndims,1]));

end
































