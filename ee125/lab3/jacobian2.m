% This program calculates the spatial Jacobian for the
%  the AdeptSix 300 robot arm.

%first calculate jacobian
%second run through all values of theta (for loop), checking where
%  det(J) = 0.  When it is (if det(J) == 0), find the joint velocities
%  which produce no end effector velocities (i.e., the nullspace of J,
%  null(J)).  Then find the the workspace velocities that can't be
%  achieved, which is the left nullspace of J, or the nullspace of the
%  transpoose of J: null(J').

clear;clc;
big = 100000;     % Checking for singularity...

w(:,1) = [0 0 1]'; q(:,1) = [0 0 0]';
w(:,2) = [0 1 0]'; q(:,2) = [150 0 0]';
w(:,3) = [0 1 0]'; q(:,3) = [150 0 260]';
w(:,4) = [1 0 0]'; q(:,4) = [410 0 320]';
w(:,5) = [0 1 0]'; q(:,5) = [410 0 320]';
w(:,6) = [0 0 1]'; q(:,6) = [410 0 320]';

% for th(1) = linspace(-2.967,2.967,1000)
% for th(2) = linspace(-0.7854,2.618,1000)
% for th(3) = linspace(-3.316,1.222,1000)
% for th(4) = linspace(-2*pi,2*pi,1000)
% for th(5) = linspace(-0.7854,3.9270,1000)
% for th(6) = linspace(-2*pi,2*pi,1000)



for i = 1:6
    v(:,i) = cross(-w(:,i),q(:,i));
    wprime(:,i) = w(:,i);
    for j = i-1:-1:1
        SE3 = screw(w(:,j),v(:,j),1);
        % obviously you have to change th = 1
        wprime(:,i) = rodrigues(w(:,j),1)*wprime(:,i);



        xi(:,i) = [v(:,i); w(:,i)];
        Jst(:,:,i) = wedge(xi(:,i));
        for j = i-1:-1:1
            SE3 = screw(w(:,j),v(:,j),1);
            Jst(:,:,i) = SE3*Jst(:,:,i)*inv(SE3);
        end
        J(:,i) = vee(Jst(:,:,i));
    end
    J
    % if cond(J) > big
    %     fprintf('Singular!\n
    % end
