function logsum = logprobs2logcumsumprobs(b)
% logprobs2logcumsumprobs  Log of the cumsum of exponentiated small numbers
%
%   USAGE:
%       logcmf = logprobs2logcumsumprobs(logpmf);
%
%   logprobs2logcumsumprobs transforms log probabilities into the log of 
%   the cumulative sum of the probabilities, i.e., 
%           LOGSUM = log(cumsum(exp(b_k))).
%   But it doesn't do the calculation so stupidly, since, if the elements
%   of b are log-probabilities, that equation would cause rounding errors.
%   The true calculation follows notes by Roweis ("Unsupervised Learning &
%   EM Algorithm," 2003).
%
%   NB: each *column* is a separate datum.
%   
%   Cf. logprobs2logsumprobs.m.

%-------------------------------------------------------------------------%
% Cribbed: 08/05/16
%   -from logprobs2logsumprobs.m
%   by JGM (but see above)
%-------------------------------------------------------------------------%

MAXEXPONENT = log(realmax);

B = MAXEXPONENT - log(size(b,1)) - max(b) - 1;
%%% I shouldn't have to subtract 1, but somehow in practice I do.  But
%%% don't worry, the calculation is correct no matter what B is; it's just
%%% that the bigger B, the better.
logsum = log(cumsum(exp(b + B))) - B;

end