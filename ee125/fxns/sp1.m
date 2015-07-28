function theta = sp1(w,u,v);
% SP1 computes the inverse kinematics Paden-Kahan subproblem 1.
%  Given the axis of rotation w, a vector u connecting a point on the axis
%  to the initial location, and a vector v connecting that point on the
%  axis to the final location (see Fig. 3.8 in Murray, Lee, Sastry), it
%  returns the value of theta, the angle of rotation.
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

uprime = u - w*w'*u;
vprime = v - w*w'*u;
theta = atan2(w'*(cross(uprime,vprime)),dot(uprime,vprime));