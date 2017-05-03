function IntegL = getInitialIntegs(Di,Si,wts,params)
% grabs integrated estimates in local space, from this dataset (Di)

%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -didn't actually change anything---but shortdata and estStatsCorePP now
%   worth with different-sized (tensor) objects shatL.
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

Do = updownDBN(Di,wts,params,'means'); % 'Nsamples');
[~, ~, shatL, ~] = estStatsCorePP(Si,params,'CoM',Do);
    
end