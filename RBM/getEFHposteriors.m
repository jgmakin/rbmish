function [p,q,r] = getEFHposteriors(R0,wts,params)
%%%%% may want to rename this
% Compute the posterior cumulants for the unisensory

%-------------------------------------------------------------------------%
% Revised: 07/08/14
%   -re-wrote from scrach to use tensor operations rather than loops.
%   -outputs are now structures with tensor fields, rather than structure
%   arrays with matrix fields
%   -incorporated posteriorEtaBars.m
%   -renamed: getSuffStats.m -> getEFHposteriors.m
% Revised: 12/17/12
%   -parallelized
% Cribbed: 12/13/12
%   from posteriorRecovery.m
%   by JGM
%-------------------------------------------------------------------------%

% init
p = getPosteriorCumulants(R0,params);

[~,R1] = updown(R0,wts,params,'propagation','Nsamples');
q = getPosteriorCumulants(R1,params);

[~,R2] = updown(R0,wts,params,'propagation','means');
r = getPosteriorCumulants(R2,params);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function p = getPosteriorCumulants(R,params)

% init
N = params.N;
C = params.C;
g = params.g;
Ndims = params.Ndims;
Nmods = params.Nmods;
gridsize = params.gridsize;
Nexamples = size(R,1);

% compute p(s|r1), p(s|r2) (etc.), p(s|r1,r2)
[cntrOfMass, ttlSpks] = GTPNsuffstats(R,params);
Info = GTPNposteriorInfo(ttlSpks,params);
unisensoryCumulants = cumulantNeutralize(cntrOfMass,Info,params);
multisensoryCumulants = gaussPosteriorization(unisensoryCumulants);

% store
p = unisensoryCumulants;
p.Xpct(:,:,end+1) = multisensoryCumulants.Xpct;
p.Info(:,:,:,end+1) = multisensoryCumulants.Info;
p.ttlSpks = ttlSpks;
p.srcs{end+1} = 'opt';

% get the g -> eta scalar
rho = ((N-1)/gridsize)^Ndims;
rhoZ = rho*((2*pi)^(Ndims/2))*sqrt(det(C));

% compute p(s|eta1bar,eta2bar,psi1,psi2)
ttlSpksBar = rhoZ*g*ones(Nexamples,Nmods);
Info = GTPNposteriorInfo(ttlSpksBar,params);
unisensoryCumulants = cumulantNeutralize(cntrOfMass,Info,params);
multisensoryCumulants = gaussPosteriorization(unisensoryCumulants);

% store
p.Xpct(:,:,end+1) = p.Xpct(:,:,end);
%%% just set to E[s|v], but you should only use this in the calculation
%%% of KL divergence that ignores the contribution of mean differences
p.Info(:,:,:,end+1) = multisensoryCumulants.Info;
p.srcs{end+1} = 'optetabar';

end
%-------------------------------------------------------------------------%
