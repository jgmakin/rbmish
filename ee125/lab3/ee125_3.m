% Lab 3: AdeptSix 300 Lab
%-------------------------------------------------------------------------%
% This program first calculate the inverse kinematics solution
%  set, [th1 th2 th3 th4 th5 th6], given the desired end effector position
%  x,y,z; and the yaw, pitch, and roll y,p,r; in the vector position =
%  [x,y,z,y,p,r].  Angles are given in rad/s.  
% It then calculates the end-effector position as a check, 
%  using the joint angles just found.
%-------------------------------------------------------------------------%
% Files: invkin.m, sp1.m, sp2.m, sp2prime.m, screw.m, ypr.m,
%  rodrigues.m, jacobian.m, wedge.m, vee.m, and ssm.m
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Revised: 11/12/03
%   -last revision for ee125
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%




clear;clc;
% gd = [-0.4317 -0.5344 0.7266 325.3; 0.5152 -0.8073 -0.2877 92.3;...
%        0.7404 0.2502 0.6239 508.6; 0 0 0 1 ]
% For testing purposes:
% position = [500 0 320 0 pi/2 0];
% position = [560.122, 0, -47.573,0,pi,0]
% position = [400.0, -300.0, -42.0,0,pi,0]

position = [-400,-350,0,90,180,250]
position(4:6) = position(4:6)/180*pi;

gst = [eye(3) [150+260; 0; 260+60+90]; 0 0 0 1];

w(:,1) = [0 0 1]'; q(:,1) = [0 0 0]';
w(:,2) = [0 1 0]'; q(:,2) = [150 0 0]';
w(:,3) = [0 1 0]'; q(:,3) = [150 0 260]';
w(:,4) = [1 0 0]'; q(:,4) = [410 0 320]';
w(:,5) = [0 1 0]'; q(:,5) = [410 0 320]';
w(:,6) = [0 0 1]'; q(:,6) = [410 0 320]';
for i = 1:6
    v(:,i) = cross(-w(:,i),q(:,i));
end

%forward kinematics check

% % For testing purposes
 p = position(1:3)';
 R = ypr(position(4),position(5),position(6));
 gd = [R p; 0 0 0 1];


theta = invkin(gd);

for j = 1:size(theta,1)
    gcheck(:,:,j) = gst;
    for i = 6:-1:1
        gcheck(:,:,j) = screw(w(:,i),v(:,i),theta(j,i))*gcheck(:,:,j);
    end
    gcheck(:,:,j)
end    

%Jacobian/singularity check
J = manjacobian(0, 0, 0, 0, pi/2, 0)
% This configuration makes axes 4 and 6 colinear.  Two colinear
%  revolute joints result in a singularity.
null(J)
% This is an orthonormal basis for the joint velocities which produce
%  no end effector velocities.         [0 0 0 sqrt(2)/2 0 -sqrt(2)/2]
null(J')
% This is an orthonormal basis for the workspace velocities that
%  can't be achieved.                  [0 0 0 1 0 0]