function [mu0p,cov0p] = priorcreator(mu,numstd,params)
% [mu0v cov0v mu0p cov0p] = priorcreator(mu,numstd,params)
% 
% Outputs the mean (mu0p) and variance (cov0p) of a distribution s.t.:
%
%   mu0p:(smin,smax)::mu:(0,1)
%   there are numstd standard deviations b/n smin and smax
%
% The smin and smax are chosen based on params.NS.  
%
% NB that Gaussian priors in 2D vis space are not very useful.


%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -changed to create prior in space of params.NS, rather than just prop
%   space.  
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


% init
sminNS = params.smin(:,strcmp(params.mods,params.NS));
smaxNS = params.smax(:,strcmp(params.mods,params.NS));
Ndims = params.Ndims;

% compute
if numstd > 0
    Qp = diag(smaxNS - sminNS)*eye(params.Ndims)/numstd;  % std matrix
    cov0p = Qp*Qp';
else
    cov0p = zeros(Ndims);
end

mu0p = scalefxn(mu',zeros(Ndims,1),ones(Ndims,1),sminNS,smaxNS)';