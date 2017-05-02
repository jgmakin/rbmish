function logsum = logprobs2logsumprobs(b)
% logprobs2logsumprobs  Log of the sum of exponentiated small numbers
%   logprobs2logsumprobs transforms log probabilities into the log of the
%   sum of the probabilities, i.e., 
%
%           LOGSUM = log(sum(exp(b_k))).
%
%   But it doesn't do the calculation so stupidly, since, if the elements
%   of b are log-probabilities, that equation would cause rounding errors.
%   The true calculation follows notes by Roweis ("Unsupervised Learning &
%   EM Algorithm," 2003).
%
% NB: each *column* is a separate datum.

%-------------------------------------------------------------------------%
% Revised: 07/13/16
%   -augmented to run on matrix inputs
% Created: 03/20/12
%   by JGM (but see above)
%-------------------------------------------------------------------------%

MAXEXPONENT = log(realmax);

B = MAXEXPONENT - log(size(b,1)) - max(b) - 1;
%%% I shouldn't have to subtract 1, but somehow in practice I do.  But
%%% don't worry, the calculation is correct no matter what B is; it's just
%%% that the bigger B, the better.
logsum = log(sum(exp(b+B))) - B;

end