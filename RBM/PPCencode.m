function R = PPCencode(X,g,inputUnitType,params)
% PPCencode   Response of a population of Gaussian-tuned neurons
%
% USAGE:
%   R = PPCencode(X,g,inputUnitType,params)
%
% Samples a matrix of size (Nexamples x Nneurons), when given a matrix of
% inputs X of size (Ndims x Nexamples) (notice this transposition), the
% type of neuron (inputUnitType), the tuning curve maxima g (Nexamples x 1)
% and the params structure.

%-------------------------------------------------------------------------%
% Revised: 02/17/15
%   -replaced call to tiledGaussianTuning.m with call to tiledTuning.m
%   -renamed from GTrespfxn.m to PPCencode.m
% Revised: 06/21/14
%   -put the guts into tiledGaussianTuning
%   -took out the sampling of gains, which now lives in
%   unifSmplAboutCenter.m, and should be called before this function
% Revised: 05/06/14
%   -vectorized inputs and outputs! (avoids parfor loop)
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% embed in smoothly tiled Gaussian tuning curves and sample
F = tiledTuning(X,g,inputUnitType,params);
R = sampler(F,inputUnitType,params);

end
