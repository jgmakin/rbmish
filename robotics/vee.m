function xi = vee(xihat)
% VEE   The 'v' operator
% VEE(xi) returns the six-dimensional twist parameterization, xi = [v; w], 
%   given a matrix xihat, a member of se(3) in homogeneous coordinates.
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -pushed one line into a separate fxn
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

omegahat = xihat(1:3,1:3);
v = xihat(1:3,4);

% may want to replace this in the future w/an "antissm" function
omega = so32rotAxis(omegahat);

xi(1:3,1) = v;
xi(4:6,1) = omega;

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function omega = so32rotAxis(omegahat)
% so32rotAxis   Transform from so(3) to a rotation axis
% so32rotAxis returns the vector omega, given the 3x3 skew-symmetric matrix
% omegahat, where omegahat*b = cross(omega,b) for any b.  See rotAxis2so3.m
% for more details.

omega(1,1) = omegahat(3,2); omega(2,1) = omegahat(1,3); omega(3,1) = omegahat(2,1);
rotAxis2so3
end
%-------------------------------------------------------------------------%