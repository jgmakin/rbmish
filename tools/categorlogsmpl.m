function ind = categorlogsmpl(logpmf,n)
% CATEGORLOGSMPL   Sampler from a categorical distribution
%   CATEGORLOGSMPL(LOGPMF,N) draws N samples from a categorical 
%   distribution with *log* probabilities LOGPMF, where exp(LOGPMF) ought 
%   to sum to one.

%-------------------------------------------------------------------------%
% Created: 03/20/12
%   by JGM
%-------------------------------------------------------------------------%

% get the log cumulative mass function
logcmf = zeros(size(logpmf));
for i = 1:length(logpmf)
    logcmf(i) = logprobs2logsumprobs(logpmf(1:i));
end

M = repmat(logcmf,1,n)>repmat(log(rand(1,n)),length(logcmf),1);
[ind, ~] = find(diff([zeros(1,n);M]));

end