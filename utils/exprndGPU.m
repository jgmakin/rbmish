function smpls = exprndGPU(mu)
% exprndGPU     Sample from an exponential distribution, using the GPU

%-------------------------------------------------------------------------%
% Created: 07/29/14
%   by JGM
%-------------------------------------------------------------------------%

mu(mu < 0) = NaN;
smpls = -mu.* log(gpuArray.rand(size(mu)));

end