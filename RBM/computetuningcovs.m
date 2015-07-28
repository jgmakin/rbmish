function tuningCov = computetuningcovs(params)
% COMPUTETUNINGCOVS     Tuning-curve covariances

%-------------------------------------------------------------------------%
% Cribbed: 02/14/12
%   -from gainerrors.m (you can use it in a lot of other places, too)
%   by JGM
%-------------------------------------------------------------------------%


tuningCov = cell(params.Nmods,1);
range = params.smax - params.smin;
C = params.C;
respLength = params.respLength;
for iMod = 1:params.Nmods
    tuningCov{iMod} = (diag(range(:,iMod)/respLength)*sqrtm(C))^2;
end

end