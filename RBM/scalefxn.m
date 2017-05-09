function Z = scalefxn(X,xmin,xmax,zmin,zmax)
% SCALEFXN  Rescales variables
%   Z = scalefxn(X,xmin,xmax,zmin,zmax)
%
%   SCALEFXN transforms x \in [xmin, xmax] into z \in [zmin, zmax] using an
%   affine rescaling.
%
%   NB: X must have size Nexamples x Ndims, where Ndims is the common
%   length of xmin, xmax, zmin, and zmax.  That is, X is in JGMSTD format.

%-------------------------------------------------------------------------%
% Revised: 08/24/16
%   -changed expected shape of inputs!!
% Revised: 07/24/14
%   -updated to use binary singleton expansion (bsxfun)
% Revised: 05/30/13
%   -made it work for both vector and matrix inputs
% Revised: 11/09/10
%   -actually complete rewritten: abolished for-loop in favor of vector
%   calculations
% Created: 07/05/10
%   by JGM
%-------------------------------------------------------------------------%

% vectorize
xmin = xmin(:)'; zmin = zmin(:)'; 
xmax = xmax(:)'; zmax = zmax(:)'; 

% lots of implicit expansion
colScaling = (zmax - zmin)./(xmax - xmin);
Z = colScaling.*(X - xmin) +  zmin;

end