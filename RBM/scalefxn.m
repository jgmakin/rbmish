function Z = scalefxn(X,xmin,xmax,zmin,zmax)
% SCALEFXN  Rescales variables
%   Z = scalefxn(X,xmin,xmax,zmin,zmax)
%
%   SCALEFXN transforms x \in [xmin, xmax] into z \in [zmin, zmax] using an
%   affine rescaling.
%
%   NB that X is expected to have the same number of *rows* as elements of
%   xmin (or xmax, zmin, zmax).

%-------------------------------------------------------------------------%
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
xmin = xmin(:); zmin = zmin(:); 
xmax = xmax(:); zmax = zmax(:); 

% opaque, fast code
rowScaling = (zmax - zmin)./(xmax - xmin);
Z = bsxfun(@plus, bsxfun(@times,rowScaling,bsxfun(@minus,X,xmin)), zmin);

end