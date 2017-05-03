function g = screw(varargin)
% screw     Screw motion from axis, pitch, and magnitude
% 
% USAGE:
%   g = screw(w,v);
%   g = screw(w,v,th);
%
% Given the axis of rotation w (size 3x1) and pitch v (size 3x1), i.e. the
% twist coordinates, returns the corresponding screw motion g.  By default,
% the magnitude of the screw rotation is assumed to be 1, but it may also
% be specified by a third argument, th.  The screw motion g is the 
% homogeneous representation of SE(3) corresponding to the twist
% coordinates.

%-------------------------------------------------------------------------%
% Revised: 11/17/16
%   -rewrote "help"
% Revised: 07/02/14
%   -changed if statement to use ~any rather than checking every entry
%   -vectorized input, tensorized output
% Revised: 08/19/10
%   -changed to varargin for backward compatibility with no th input
%   version ("twist.m")
% Revised: 06/09/10
%   -copied exactly from twist2.m (I think it ought to have been called
%   "screw" all along)
% Revised: 06/09/10
%   -commented
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

% collect inputs
omega = varargin{1};
v = varargin{2};
if length(varargin) < 3, th = 1; else th = varargin{3}; end

% Ns
Nangles = length(th);
Ndims = length(omega);

% rotation matrices via Rodrigues's formula
R = rodrigues(omega, th);

% now compute SE(3) matrices
th = reshape(th(:),[1,1,Nangles]);
if ~any(omega)                                  % prismatic joints
    p = v.*th;
else                                            % revolute joints
    % p = ((eye(3) - R)*omegahat + omega*omega'*th)*v;
    omegahat = rotAxis2so3(omega);
    p = eye(3)*omegahat*v + ...
        (-tensorOp(R,repmat(omegahat*v,[1,1,Nangles])) + (omega*omega'*v).*th);
    %%%% may want to rethink this ^, using implicit expansion 
end
% g = [R, p; 0 0 0 1];
g = cat(1,cat(2,R,p),cat(2,zeros(1,Ndims,Nangles),ones(1,1,Nangles)));


end






