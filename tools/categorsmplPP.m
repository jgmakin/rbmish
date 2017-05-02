function ind = categorsmplPP(pmf,n)
% CATEGORSMPLPP   Sampler from a categorical distribution
%   CATEGORSMPLPP(PMF,N) draws N samples from a categorical distribution 
%   with probabilities PMF, where PMF ought to sum to one.
%
%   This version uses the parallel-processing toolbox!  It's better for
%   large n (~40000).

%-------------------------------------------------------------------------%
% Created: 03/20/12
%   by JGM
%-------------------------------------------------------------------------%

cmf = cumsum(pmf);

ind = zeros(n,1,'like',pmf);
parfor i = 1:n, 
    ind(i) = find(cmf>rand,1);
end

end
