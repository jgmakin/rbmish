function [M, omega] = rotation(x);
%ROTATION calculates the composite rotation matrix, M, and the normalized 
% equivalent axis, omega, for the following rotations around a fixed, 
% spatial frame for small values of x:
%
% (a) Rotate around the x-axis by x radians
% (b) Rotate around the y-axis by x radians
% (c) Rotate around the x-axis by ?x radians
% (d) Rotate around the y-axis by ?x radians

XPOS = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
YPOS = [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)];
XNEG = [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)];
YNEG = [cos(x) 0 -sin(x); 0 1 0; sin(x) 0 cos(x)];
M = YNEG*XNEG*YPOS*XPOS;

theta = acos((trace(M)-1)/2);
omega = 1/(2*sin(theta))*[M(3,2)-M(2,3); M(1,3)-M(3,1); M(2,1)-M(1,2)];
omega = omega/sum(omega.^2);