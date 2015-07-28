function varargout = longdata(varargin)
% Go from the standard tensor to one with all examples in the first
% dimension.  The order of the (cases? batches?) is preserved.
%
% USAGE:
%   Slong = longdata(Sshort)
%   [Slong,Ylong] = longdata(Sshort,Yshort)

%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed to work with tensors with more than three "dimensions"---as
%   long as the first and last are the number of cases and batchsize!
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


for i = 1:nargin
    shrtform = varargin{i};
    
    shrtdims = size(shrtform);
    longdims = [shrtdims(2:end-1) shrtdims(1)*shrtdims(end)];
    cmptform = reshape(shiftdim(shrtform,1),longdims);
    longform = shiftdim(cmptform,length(longdims)-1);
    
    varargout{i} = longform;
end

end
