function g = screw(varargin)
%SCREW takes as input two 3×1 vectors w and v as twist coordinates, and
%   calculates the homogeneous representation of SE(3) corresponding to
%   these coordinates.
%-------------------------------------------------------------------------%
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
w = varargin{1};
v = varargin{2};
if length(varargin) < 3, th = 1; else th = varargin{3}; end


Nangles = length(th);
Ndims = length(w);

R = rodrigues(w, th);
if ~any(w)                                      % prismatic joints
    p = v*th(:)';
else                                            % revolute joints
    what = ssm(w);
    % p = ((eye(3) - R)*what + w*w'*th)*v;          
    p = repmat(what*v,[1,Nangles]) -...
        reshape(reshape(R,[Ndims*Nangles,Ndims])*what*v,[Ndims,Nangles]) +...
        w*w'*v*th(:)';
end
% g = [R, p; 0 0 0 1];
g = cat(3,cat(1,R,zeros(1,Nangles,Ndims)),cat(1,p,ones(1,Nangles)));

end






