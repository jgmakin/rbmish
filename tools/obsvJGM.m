function OB = obsvJGM(A,C,varargin)
% obsvJGM   Observability matrix
%
% USAGE:
%
%   OB = obsvJGM(A,C);          % observability matrix
%   OB = obsvJGM(A,C,100)       % extended (to 100) observability matrix
%
% This replaces Matlab's obsv.m---so you don't have to use the control
% toolbox.  It also allows construction of the "extended" observability
% matrix, i.e. up to the Nth power (rather than n, the dimension of the
% state), via a variable input argument.

%-------------------------------------------------------------------------%
% Created: 10/16/14
%   by JGM
%-------------------------------------------------------------------------%


if isempty(varargin), N = size(A,1); else N = varargin{1}; end

OB = C;
for n = 1:N-1, OB = cat(1,C,OB*A); end
    
    
end