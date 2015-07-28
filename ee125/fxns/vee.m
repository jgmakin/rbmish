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

what = xihat(1:3,1:3);
v = xihat(1:3,4);

% may want to replace this in the future w/an "antissm" function
w = antissm(what);

xi(1:3,1) = v;
xi(4:6,1) = w;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function w = antissm(what);
% ANTISSM   Transforming a skew-symmetric matrix to a vector
% ANTISSM returns the vector w, given the 3x3 skew-symmetric matrix what,
% where what*b = cross(w,b) for any b.  See ssm.m for more details.
%-------------------------------------------------------------------------%
% Created: 06/09/10
%   by JGM
%-------------------------------------------------------------------------%

w(1,1) = what(3,2); w(2,1) = what(1,3); w(3,1) = what(2,1);

end