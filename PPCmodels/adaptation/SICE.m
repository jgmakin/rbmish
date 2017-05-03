function [SINSCerrMu,SINSCondCov] = SICE(s,params,varargin)
% SICE  single-input conditional error stats


%-------------------------------------------------------------------------%
% Revised: 12/11/13
%   -changed the order in which variables are stored in the output
%   arguments to the "new standard" (see setParams.m)
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% init
Ndims = params.Ndims;
Nmods = length(params.mods); % 2;
n = nargin - 2;

% malloc
SINSCondCov = zeros(Ndims,Ndims,Nmods,n); 
SINSCerrMu = zeros(Ndims,Nmods,n);

% convert to "neutral" space
for iMod = 1:Nmods
    
    J = ntrlJacobian(s,iMod,params);
    for iArg = 1:n
        SINSCondCov(:,:,iMod,iArg) = J*varargin{iArg}{iMod}.cov*J';
        SINSCerrMu(:,iMod,iArg) = J*varargin{iArg}{iMod}.mu;
    end
end

end