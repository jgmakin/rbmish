function d = covdist(M,P)
% COVDIST   Distance between two positive semidefinite matrices
%   COVDIST computes the distance b/n two covariance matrices M and P 
%   according to the appropriate distance measure, 
%
%           norm(log(spectrum(M\P))).
%
%   (See Mitteroecker 2009.)

%-------------------------------------------------------------------------%
% Created: 01/30/12
%   by JGM
%-------------------------------------------------------------------------%

d = norm(log(eig(M\P)));

end