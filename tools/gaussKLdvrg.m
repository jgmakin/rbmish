function KL = gaussKLdvrg(mu0,mu1,Prcn0,Prcn1)
% gaussKLdvrg   KL divergence between two normal distributions
%
% USAGE:
% 
%   KL = gaussKLdvrg(mu0,mu1,Prcn0,Prcn1)
% 
% The KL divergence b/n two normal distributions has a simple formula---see
% the wikipedia entry on the MND, e.g., or rederive it!
%
% NB that this is the divergence from N0 to N1: KL(N0,N1)
%
% Also NB that the third and fourth args are PRECISIONS!!

%-------------------------------------------------------------------------%
% Created: 10/22/12
%   by JGM
%-------------------------------------------------------------------------%


if size(mu0,1)~=length(mu0), mu0 = mu0'; end
if size(mu1,1)~=length(mu1), mu1 = mu1'; end

S = Prcn1/Prcn0;
m = (mu1-mu0);
KL = (trace(S) + m'*Prcn1*m - log(det(S)) - length(m))/2;

end

%%%%%%%
% det(S)    = det(P1*inv(P0))
%           = det(P1)*det(inv(P0))
%           = det(P1)/det(P0)
%
% tr(S)     = tr(P1*inv(P0))
%           = tr(inv(P0)*P1)
%           = 
%%%%%%%
