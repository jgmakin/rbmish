function [Y,SigmaY] = PPC2FltrCmlnts(PPCs,S,Z,params)
% PPC2FltrCmlnts    Probabilistic population codes --> filter cumulants
%
% USAGE:
%
%       [Y,SigmaY] = PPC2FltrCmlnts(PPCs,S,Z,params)
%
% PPC2FltrCmlnts takes PPCs from the recurrent EFH; the encoded stimuli, S;
% the state, Z; and the parameters structure; and returns the corresponding
% terms in the classic linear dynamical system, Y and SigmaY.  These are
% *the emission* and the *the emission covariance*, resp.
%
% NB: PPCs, S, and Z should be in shortdata (standard) form.

%-------------------------------------------------------------------------%
% Revised: 08/29/14
%   -fixed a bug, making Z an argument
% Created: 08/26/14
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Ncases,Ndims,Nmods,T] = size(S);

% get the posterior (assuming flat prior) cumulants
[cntrsOfMass, ttlSpks] = GTPNsuffstats(longdata(PPCs),params);
Y = shortdata(Ncases,4,cntrsOfMass);
Cvrn = shortdata(Ncases,5,GTPNposteriorCvrn(ttlSpks,params));

% Assume no correlations b/n emissions from different populations
SigmaY = zeros(Ncases,Ndims*Nmods,Ndims*Nmods,T,'like',S);
for iMod = 1:Nmods
    inds = 1+(iMod-1)*Ndims:(iMod*Ndims);
    SigmaY(:,inds,inds,:) = Cvrn(:,:,:,iMod,:);
end

% if necessary, correct the emissions in the case of a torus
if strcmp(params.dynamics.walls,'wrapping')
    if isfield(params.dynamics,'H') %%% not ideal
        C = blkdiag(params.dynamics.C,params.dynamics.H);
    else
        C = params.dynamics.C;
    end
    Y = getBestTorusY(Y,S,Z,C,params);
else
    Y = reshape(Y,[Ncases,Ndims*Nmods,T]);
end

end










