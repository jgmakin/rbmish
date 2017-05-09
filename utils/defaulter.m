function val = defaulter(name,defaultval,varargin)
% defaulter     Assign default or variable-input value
%
% USAGE:
%   val = defaulter(name,defaultval,varargin{:})
% 
% Often you want to have a default value for some variable, but allow it to
% take on some new value based on a variable input argument.  Matlab's
% standard way of dealing with this in its own functions is to have the
% user specify the name, followed by the value, of the variable, in
% consecutive arguments to a function, usually at the end.  This function
% implements such a scheme when placed inside some function
%
%   foo(a,b,...,varargin)
%
% according to the usage above. 

%-------------------------------------------------------------------------%
% Created: 09/29/16
%   by JGM
%-------------------------------------------------------------------------%


iVal = find(strcmp(varargin,name));
if ~isempty(iVal), val = varargin{iVal+1}; else val=defaultval; end


end