function ind = categorlogsmplPP(logpmf,n)
% CATEGORLOGSMPL   Sampler from a categorical distribution
%   CATEGORLOGSMPL(LOGPMF,N) draws N samples from a categorical 
%   distribution with *log* probabilities LOGPMF, where exp(LOGPMF) ought 
%   to sum to one.
%
%   This version uses the parallel-processing toolbox!  It's better for
%   large n (~40000).

%-------------------------------------------------------------------------%
% Revised: 08/05/16
%   -replaced parfor loop with loopless function logprobs2logcumsumprobs.m
%       to convert logpmf into logcmf.
% Created: 03/20/12
%   by JGM
%-------------------------------------------------------------------------%

% get the log cumulative mass function
logcmf = logprobs2logcumsumprobs(logpmf);

% sample from it
parfor i = 1:n
    ind(i) = find(logcmf>log(rand),1);
end

end