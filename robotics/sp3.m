function theta = sp3(w,p,q,r,del)
% SP3 computes the inverse kinematics Paden-Kahan subproblem 3.
%  This function calculates the angle of rotation from a point p to a point
%  a distance delta away from the final point q, given p, q, delta, and the
%  axis of rotation w (see Fig. 3.10 in Murray, Lee, Sastry).
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

u = p - r;
v = q - r;
uprime = u - w*w'*u;
vprime = v - v*v'*v;
delprimesq = del^2 - abs(w'*(p-q))^2;   %????
th0 = atan2(w'*cross(uprime,vprime),dot(uprime,vprime));
theta = th0 + acos(((norm(uprime))^2 +...
    (norm(vprime))^2 - delprimesq)/(2*norm(uprime)*norm(vprime)));
theta(2,1) = th0 - acos(((norm(uprime))^2 +...
    (norm(vprime))^2 - delprimesq)/(2*norm(uprime)*norm(vprime)));