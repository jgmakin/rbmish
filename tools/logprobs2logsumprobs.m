function logsum = logprobs2logsumprobs(b)
% LOGPROBS2LOGSUMPROBS  Log of the sum of exponentiated small numbers
%   LOGPROBS2LOGSUMPROBS transforms log probabilities into the log of the
%   sum of the probabilities, i.e., 
%           LOGSUM = log(sum(exp(B_k))).
%   But it doesn't do the calculation so stupidly, since, if the elements
%   of B are log-probabilities, that equation would cause rounding errors.
%   The true calculation follows notes by Roweis ("Unsupervised Learning &
%   EM Algorithm," 2003).

%-------------------------------------------------------------------------%
% Created: 03/20/12
%   by JGM (but see above)
%-------------------------------------------------------------------------%

MAXEXPONENT = log(realmax);

B = MAXEXPONENT - log(length(b)) - max(b) - 1;
%%% I shouldn't have to subtract 1, but somehow in practice I do
logsum = log(sum(exp(b+B))) - B;

end