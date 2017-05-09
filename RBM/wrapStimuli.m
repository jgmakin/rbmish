function [S,srange] = wrapStimuli(S,smin,smax,N)
% wrapStimuli   Wrap stimuli onto an N-torus
%
% USAGE:
%   S = wrapStimuli(S,smin,smax,N);
%
%   [S,srange] = wrapStimuli(S,smin,smax,N);


%-------------------------------------------------------------------------%
% Cribbed: 01/05/17
%   -from getDataPPC
%   -by JGM
%-------------------------------------------------------------------------%

% Ns
[~,Ndims,Nmods] = size(S);

% now wrap the stimuli into encoding space
smin = reshape(smin,[1,Ndims,Nmods]);       % reshape for implicit expans.
smax = reshape(smax,[1,Ndims,Nmods]);       %
srange = N/(N-1)*(smax - smin);             % confusing, to be sure
S = mod(S-smin,srange) + smin;

end