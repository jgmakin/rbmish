function xihat = wedge(xi)
% WEDGE   The 'wedge' operator
% WEDGE(xi) returns the matrix xihat, a member of se(3), in homogeneous
%  coordinates, given its six-dimensional parameterization, xi = [v; w].

%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

w = [xi(4); xi(5); xi(6)];
v = [xi(1); xi(2); xi(3)]; 
xihat = [rotAxis2so3(w) v; 0 0 0 0];