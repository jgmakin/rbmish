function ind = categorlogsmpl(logpmf,n)
% categorlogsmpl   Sampler from a categorical distribution
%   categorlogsmpl(logpmf,n) draws n samples from a categorical 
%   distribution with *log* probabilities logpmf, where exp(logpmf) ought 
%   to sum to one.

%-------------------------------------------------------------------------%
% Revised: 08/05/16
%   -replaced loop with loopless function logprobs2logcumsumprobs.m to 
%       convert logpmf into logcmf.
% Created: 03/20/12
%   by JGM
%-------------------------------------------------------------------------%

% get the log cumulative mass function
logcmf = logprobs2logcumsumprobs(logpmf);

% sample from it 
M = logcmf > log(rand(1,n));
[ind, ~] = find(diff([zeros(1,n);M]));

end
