function J = manjacobian2(th1,th2,th3,th4,th5,th6)
% MANJACOBIAN   Manipulator spatial jacobian
%   This program calculates the spatial Jacobian for the AdeptSix 300 robot 
%   arm, given a joint configuration
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

w(:,1) = [0 0 1]'; q(:,1) = [0 0 0]';
w(:,2) = [0 1 0]'; q(:,2) = [150 0 0]';
w(:,3) = [0 1 0]'; q(:,3) = [150 0 260]';
w(:,4) = [1 0 0]'; q(:,4) = [410 0 320]';
w(:,5) = [0 1 0]'; q(:,5) = [410 0 320]';
w(:,6) = [0 0 1]'; q(:,6) = [410 0 320]';

th = [th1 th2 th3 th4 th5 th6];

for i = 1:6
    v(:,i) = cross(-w(:,i),q(:,i));
    xi(:,i) = [v(:,i); w(:,i)];
    Jst(:,:,i) = wedge(xi(:,i));   
    for j = i-1:-1:1
        SE3 = screw(w(:,j),v(:,j),th(j));
        Jst(:,:,i) = SE3*Jst(:,:,i)*inv(SE3);  
    end
    J(:,i) = vee(Jst(:,:,i));
end
J;