function [MSE, Zbar] = testEFHNextFrameError(R,X,Q,wts,params)
% getNextFrameError
%
% USAGE:
%   [R,X,Q] = params.getTestError(dataclass);
%   MSE = getNextFrameError(R,X,Q,wts,params);
%
%   [MSE,Zbar] = getNextFrameError(R,X,Q,wts,params);
%
% NB: the inputs are in longdata format, but Zbar is shortdata!


%-------------------------------------------------------------------------%
% Revised: 12/07/16
%   -rewrote: separated out the "next-frame-error" functionality (here)
%   from the trajectory-generating (confabulating) functionality (moved to
%   rEHFprobe.m)
% Revised: 12/05/16
%   -added functions for MoCap data
% Revised: 01/20/16
%   -added function for writing animated gifs
% Created: 08/10/15
%   by JGM
%-------------------------------------------------------------------------%


% where to start measuring errors
t0 = 1; % hard-coded
Nexamples = size(R,1);
Ntraj = floor(Nexamples/Q.T);

% compute hidden states ("filter")
[~,Zbar] = updownRDBN(R,wts,params,Q.T);
Zbar = shortdata(Ntraj,3,Zbar);

% now generate the set of "next frames"
Rhat = forwardGenerate(Zbar,t0,0,'inferredhidmeans',wts,params);
% %%% (1) using Gibbs *meaning* is slightly better....

% mean square error
R = shortdata(Ntraj,3,R);
SE = (Rhat(:,:,(t0+1):end) - R(:,:,(t0+1):end)).^2;
MSE = mean(SE(:));

end