function trajs = EFHdata2LDSdata(Ntraj,getLDStensor)
% EFHdata2LDSdata   Convert data for EFH.m to data for EM4LDS.m
%
% USAGE:
%    trajs = EFHdata2LDSdata(Ntraj,T,dataclass,params);

%-------------------------------------------------------------------------%
% Revised: 02/04/16
%   -edited to just set Shat = R for non-PPC data
%   -renamed PPCdata2LDSdata -> EFHdata2LDSdata
%   -moved from /PPCmodels to /RBM
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


% generate data
[Ytensor,Xtensor] = getLDStensor();
T = floor(size(Ytensor,1)/Ntraj);

% now change the shape
Y       = shortdata(Ntraj,3,Ytensor(:,:));
Ycell   = mat2cell(permute(Y,[2,3,1]),size(Y,2),T,ones(1,Ntraj));
[trajs(1:Ntraj).Y]  = Ycell{:};

Nstates = size(Xtensor,2);
Z       = shortdata(Ntraj,3,Xtensor);
Zcell 	= mat2cell(permute(Z,[2,3,1]),Nstates,T,ones(1,Ntraj));
[trajs(1:Ntraj).Z]  = Zcell{:};

% SigmaYX? or InfoYX?

end