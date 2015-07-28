function ind = categorsmpl(pmf,n)
% CATEGORSMPL   Sampler from a categorical distribution
%   CATEGORSMPL(PMF,N) draws N samples from a categorical distribution 
%   with probabilities PMF, where PMF ought to sum to one.

%-------------------------------------------------------------------------%
% Created: 03/16/12
%   by JGM
%-------------------------------------------------------------------------%

cmf = cumsum(pmf);

M = repmat(cmf,1,n)>repmat(rand(1,n),length(cmf),1);
[ind, ~] = find(diff([zeros(1,n);M]));

end