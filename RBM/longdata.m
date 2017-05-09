function varargout = longdata(varargin)
% longdata      Transform the "standard" EFH-training tensor into a matrix
% 
% USAGES:
%   Smat = longdata(Stensor)
%   [Smat,Ymat] = longdata(Stensor,Rtensor)
%
% Given a tensor Stensor of size (Nrows x Ndims x ... x Nfinal), where the
% omitted ways are not limited in number but must each have dimension 1, 
% longdata returns a matrix of size (Nfinal*Nrows x Ndims).  Thus longdata
% undoes shortdata.m.
% 
% NB! Consecutive examples in the *final* way of Stensor are preserved as
% consecutive in the rows of Smat.

%-------------------------------------------------------------------------%
% Revised: 12/07/16
%   -rewrote help
% Revised: 12/10/13
%   -changed to work with tensors with more than three "dimensions"---as
%   long as the first and last are the number of cases and batchsize!
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


for i = 1:nargin
    shrtform = varargin{i};
    Nways = ndims(shrtform);
    shrtdims = size(shrtform);
    finaldim = size(shrtform,max(Nways,3));
    longdims = [shrtdims(2:max(2,Nways-1)) shrtdims(1)*finaldim];
    cmptform = reshape(shiftdim(shrtform,1),longdims);
    longform = shiftdim(cmptform,length(longdims)-1);
    
    varargout{i} = longform;
end


end
