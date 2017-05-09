function postStats = PPCinputStats(r,tuningCov,bias)
% PPCINPUTSTATS Compute the mean and covariance of the PPC posterior
%   PPCINPUTSTATS takes an input data vector r of firing rates, breaks
%   it into nummodes populations (retrieved from the size of the bias
%   matrix), and under the assumption of smooth, dense tiling of
%   Gaussian-tuned, Poisson-noise (GTPN) neurons, computes the mean and
%   variance of the center of mass.  See integnotes.tex for the equations.

%-------------------------------------------------------------------------%
% Created: 08/12/11
%   by JGM
%-------------------------------------------------------------------------%

% init params
Nmods = size(bias,2);

% malloc
postStats = cell(Nmods,1);

eta = sum(reshape(r,length(r)/Nmods,Nmods));
for iMod = 1:Nmods
    postStats{iMod}.mu = bias(:,iMod);
    postStats{iMod}.eta = eta(iMod);
    if eta(iMod)==0
        postStats{iMod}.cov = Inf;
        % fprintf('warning: no spikes on this trial -- jgm\n');
    else 
        postStats{iMod}.cov = tuningCov{iMod}/eta(iMod);
    end
end

end