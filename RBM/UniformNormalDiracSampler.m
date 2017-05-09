function X = UniformNormalDiracSampler(xpctX,covX,Nsamples,xmin,xmax,mrgn)
% UniformNormalDiracSampler     Just like it says
%
% USAGES:
%   % generic
%   x = UniformNormalDiracSampler(xpctX,covX,Nsamples,xmin,xmax,mrgn);
%
%   % uniform
%   x = UniformNormalDiracSampler(zeros(Ndims,1),Inf,M,xmin,xmax,mrgn);
%
%   % normal
%   x = UniformNormalDiracSampler(xpctX,covX,Nsamples,[],[],[]);
%
%   % Dirac
%   x = UniformNormalDiracSampler([],0,Nsamples,[],[],[]);
% 
% Often you want to sample from a multivariate normal distribution, but
% also to accommodate as special cases sampling from a uniform distribution
% (if the covariance is infinite) or a Dirac delta (if the covariance is
% zero).  This function incorporates those three sampling functions.
%
% Notice that not all arguments are necessary for the different cases (see
% USAGES above).

%-------------------------------------------------------------------------%
% Created: 12/30/16
%   by JGM
%-------------------------------------------------------------------------%

Ndims = length(xpctX);
if sum(covX(:)) ==  Inf     % sample from a uniform distribution
    X = scalefxn(rand(Nsamples,Ndims,'like',xpctX),...
        zeros(size(xmin),'like',xpctX),ones(size(xmin),'like',xpctX),...
        xmin+mrgn,xmax-mrgn);
elseif sum(covX(:)) ==  0 % sample from a Dirac delta
    X = repmat(xpctX(:)',[Nsamples,1]);
else                        % sample from a normal distribution
    X = xpctX' + randn(Nsamples,Ndims)*chol(covX);
end

end
