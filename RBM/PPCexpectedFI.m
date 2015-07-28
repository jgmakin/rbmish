function FI = PPCexpectedFI(range,g,params)
% PPCEXPECTEDFI     Expected Fisher information, Poisson population
%   PPCEXPECTEDFI computes the Fisher information matrix FI for a
%   population of Gaussian-tuned, Poisson neurons, encoding a parameter in
%   RANGE and with additional parameters specified by PARAMS.
%-------------------------------------------------------------------------%
% Revised: 09/27/11
%   -simplified the calculation of rho*Z on the realization that it's
%       actually independent of range and response length.
% Cribbed: 11/30/10
%   from PPCinfo.m
%   by JGM
%-------------------------------------------------------------------------%

% init params
Ndims = params.Ndims;                                   % encoded vars/neuron
N = params.N;
C = params.C;
respLength = params.respLength;
gridsize0 = params.gridsize;
r = range(:,2) - range(:,1);

% rescale (if necessary)
tuningCov = (diag(r(:)/respLength)*sqrtm(C))^2;

% actually, 
rho = ((N-1)/gridsize0)^Ndims;
rhoZ = rho*((2*pi)^(Ndims/2))*sqrt(det(C));

% fisher info
FI = g*rhoZ*inv(tuningCov);

end

% D&A include a factor Ndims (dimensionality of the data), but this appears to
% be a typo: the original Kechen Zhang paper has a factor 1/Ndims.  But I have
% neither!  Check your derivation.
%
% Also note that for p=2, your favorite case, the covariance terms cancel, 
% so that Fisher Info is a fxn of gain and density alone.
%
% And finally: the density term rho is (N-1)/gridsize rather than
% N/gridsize b/c the extreme neurons land *on* the boundary (as opposed to
% deltaS/2 away from it).
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% this was your old way of calculating it; but it turned out that a lot of
% terms cancel, and rho*Z is independent of the range over which your
% encoding (and so independent of modality)

% gridsize = scalefxn(gridsize0*ones(Ndims,1),zeros(Ndims,1),respLength*ones(Ndims,1),...
%     zeros(Ndims,1),r);
% 
% % compute
% measure = prod(gridsize(:));                    % area of coverage
% rho = N^Ndims/measure;                              % neuron tiling density
% Z = (2*pi)^(Ndims/2)*sqrt(det(tuningCov));
% 
% % fisher info
% FI = g*Z*rho*inv(tuningCov);
%-------------------------------------------------------------------------%


