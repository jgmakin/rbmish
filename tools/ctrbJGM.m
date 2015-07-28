function CO = ctrbJGM(A,B,varargin)
% ctrbJGM   Controllability matrix
%
% USAGE:
%
%   CO = ctrbJGM(A,B);          % controllability matrix
%   CO = ctrbJGM(A,B,100)       % extended (to 100) controllability matrix
%
% This replaces Matlab's ctrb.m---so you don't have to use the control
% toolbox.  It also allows construction of the "extended" controllability
% matrix, i.e. up to the Nth power (rather than n, the dimension of the
% state), via a variable input argument.

%-------------------------------------------------------------------------%
% Created: 10/16/14
%   by JGM
%-------------------------------------------------------------------------%


if isempty(varargin), N = size(A,1); else N = varargin{1}; end

CO = B;
for n = 1:N-1, CO = cat(2,B,A*CO); end
    
    
end