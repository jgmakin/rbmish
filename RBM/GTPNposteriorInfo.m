function Info = GTPNposteriorInfo(ttlSpks,params)
% GTPNposteriorInfo     GTPN posterior information matrices
% USAGE:
%   Info = GTPNposteriorInfo(ttlSpks,params)
%
% Usually, ttlSpks will be provided by GTPNsuffstats.m.  NB the dimensions:
%
%   ttlSpks: (Nexamples x Nmods)
%
%   Info: (Nexamples x Ndims x Ndims x Nmods)

%-------------------------------------------------------------------------%
% Revised: 08/07/14
%   -amended for GPU data
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
Info = nan(Nexamples,Ndims,Ndims,Nmods,'like',ttlSpks);
        
% loop through modalities
for iMod = 1:Nmods
    tuningInfo = inv(tuningCov{iMod});
    Info(:,:,:,iMod) = shiftdim(reshape(...
        tuningInfo(:)*ttlSpks(:,iMod)',[Ndims,Ndims,Nexamples]),2);
end


end
