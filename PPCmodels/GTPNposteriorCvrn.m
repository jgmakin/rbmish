function Cvrn = GTPNposteriorCvrn(ttlSpks,params)
% GTPNposteriorCvrn     GTPN posterior covariance matrices
% USAGE:
%   Info = GTPNposteriorCvrn(ttlSpks,params)
%
% Usually, ttlSpks will be provided by GTPNsuffstats.m.  NB the dimensions:
%
%   ttlSpks: (Nexamples x Nmods)
%
%   Info: (Nexamples x Ndims x Ndims x Nmods)

%-------------------------------------------------------------------------%
% Revised: 07/09/14
%   -transposed ttlSpks
% Created: 07/08/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
[Nexamples,Nmods] = size(ttlSpks);
Ndims = params.Ndims;

% tuning-curve widths
tuningCov = computetuningcovs(params);

% malloc
Cvrn = NaN(Nexamples,Ndims,Ndims,Nmods,'like',ttlSpks);

% loop through modalities
for iMod = 1:Nmods
    Cvrn(:,:,:,iMod) = shiftdim(reshape(...
        tuningCov{iMod}(:)*(1./ttlSpks(:,iMod))',[Ndims,Ndims,Nexamples]),2);
end

end
