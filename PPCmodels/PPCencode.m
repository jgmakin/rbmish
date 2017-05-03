function R = PPCencode(S,g,smin,smax,dstrbs,nums,params)
% PPCencode   Response of a population of neurons to stimuli
%
% USAGE:
%   R = PPCencode(X,g,smin,smax,inputUnitType,numsUnits,params)
%
% PPCencode generates responses R (Nexamples x Nunits) to stimuli S 
% (Nexamples x Ndims), with tuning curves determined by their maxima, g
% (Nexamples x 1), the type of neuron (inputDstrb), and the parameters 
% structure (params).  The (theoretical) range of the stimulus (smin and 
% smax, with length Ndims) must also be supplied.

%-------------------------------------------------------------------------%
% Revised: 08/24/16
%   -added arguments because the rescaling now happens in tiledTuning.m
%   rather than generateData.m
%   -changed the format of input S to Nexamples x Ndims from the transpose;
%   this respects the new version of scalefxn.m
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

% embed in tuning curves and sample
F = tiledTuning(S,g,smin,smax,dstrbs,nums,params);
R = sampleT(F,dstrbs,nums,params);

end
