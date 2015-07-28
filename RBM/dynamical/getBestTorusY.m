function bestY = getBestTorusY(Y,S,Z,C,params)
% Complicated!  The point is to move the decoded Y from the toroidal
% workspace onto the true (nonwrapping) trajectory, which is in Z.
%
% This first requires computing the minimal error between the decoded (on
% the torus) Y and the true (wrapped) stimulus, S.  Then this error is
% added back into the true *non*wrapped stimulus, C*Z.

%-------------------------------------------------------------------------%
% Revised: 08/29/14
%   -now runs on both mods at once
%   -returns Y as a 3-tensor (Ncases x Ndims*Nmods x T)
% Revised: 07/07/14
%   -replaced loop calculation of CZ with tensor calculation
% Created: ??/??/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
[Ncases,Ndims,Nmods,T] = size(S);
N = params.N;

% stimulus range on a torus
smin = params.smin;
smax = params.smax;
srange = reshape(N/(N-1)*(smax - smin),[1,Ndims,Nmods]);

% compute minimal errors
eAct = Y - S;
eBck = bsxfun(@minus,eAct,srange);
eFwd = bsxfun(@plus,eAct,srange);
e = cat(5,eBck,eAct,eFwd);
[~,minInds] = min(abs(e),[],5);
[indDim1,indDim2,indDim3,indDim4] = ndgrid(1:Ncases,1:Ndims,1:Nmods,1:T);
eTrue = e(sub2ind(size(e),indDim1,indDim2,indDim3,indDim4,minInds));

% "correct" Y
eTrue = reshape(eTrue,[Ncases,Ndims*Nmods,T]);
bestY = shortdata(Ncases,3,longdata(Z)*C') + eTrue; % bestY = CZ + eTrue


end

