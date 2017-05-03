function alllattices = ndlattice(Ndims,xlv)
% ndlattice     Create an N-dimensional lattice
%
% USAGE:
%   lattice = ndlattice(Ndims,xlv)
%
% Very similar to the Mathworks' ndgrid, for which this function is a kind
% of wrapper, ndlattice returns a grid or lattice over N dimensions.  The 
% differences are that this fxn:
%   (1) uses the same lattice vector, xlv, for all N dimensions; and
%   (2) decides how many dimensions are in the lattice based *not on the
%       number of outputs, but on an argument Ndims.  
% Thus, if M = length(xlv), the output is not a set of arrays of size 
% (M x M x M x ...), but a single matrix of size M^Ndims x Ndims, each
% column of which corresponds to vectorized versions of the arrays returned
% by ndgrid.
% 
% The usefulness of this function is that it will work as is for any number
% of dimensions Ndims.

%-------------------------------------------------------------------------%
% Cribbed 01/??/17
%   -from tiledTuning.m
%   by JGM
%-------------------------------------------------------------------------%

% Ns
M = length(xlv);

lattices1Dcell = mat2cell(xlv(ones(Ndims,1,'like',xlv),:),...
    ones(Ndims,1),M);
singlelattice = ndgrid(lattices1Dcell{:});
alllattices = NaN(M^Ndims,Ndims,'like',xlv);
for iDim = 1:Ndims
    alllattices(:,iDim) = vect(shiftdim(singlelattice,iDim-1));
end


end