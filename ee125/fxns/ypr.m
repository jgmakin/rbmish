function M = ypr(y,p,r);
%YPR calculates the composite rotation matrix, M, given the yaw (y), 
% pitch (p), and roll (r), all in rad/s.  The corresponding fixed
% angle rotation matrix M is Rz(y)*Ry(p)*Rz(r).


MY = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
MP = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
MR = [cos(r) -sin(r) 0; sin(r) cos(r) 0; 0 0 1];

M = MY*MP*MR;