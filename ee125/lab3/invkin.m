function theta = invkin(gd);
% Inverse Kinematics Solution for AdeptSix 300 Robot INVKIN takes a desired
% configuration gd (a member of Special Euclidean Group 3), and returns the
% required joint angles.
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Revised: 11/18/03
%   -final revision for ee125
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%


gst = [eye(3) [150+260; 0; 260+60+90]; 0 0 0 1];
jlimhigh = [2.967, 2.618, 1.222, 1000, 3.9270, 1000];
jlimlow = [-2.967, -0.7854, -3.316, -1000, -0.7854, -1000];
index = 1;               % For use in joint limit test at the end

% Here are all the axes of rotation and their associated vectors:
w(:,1) = [0 0 1]'; q(:,1) = [0 0 0]';
w(:,2) = [0 1 0]'; q(:,2) = [150 0 0]';
w(:,3) = [0 1 0]'; q(:,3) = [150 0 260]';
w(:,4) = [1 0 0]'; q(:,4) = [410 0 320]';
w(:,5) = [0 1 0]'; q(:,5) = [410 0 320]';
w(:,6) = [0 0 1]'; q(:,6) = [410 0 320]';
for i = 1:6
    v(:,i) = cross(-w(:,i),q(:,i));
end

% ANGLE 1
% The ratio of the x and y coordinates of the point of intersection
%  of axes 4, 5, and 6 is a function only of the first joint angle.
%  Specifically, theta1 = atan(py,px).  To get this point p, we have to
%  work backwards from the end-effector position.  We do this by subtrating
%  from it the quantity 90*zt, where zt is the tool-frame z axis given in
%  spatial coordinates and 90 is the length from the point of intersection
%  to the tool frame along its z-axis.

p = gd(1:3,4) - 90*gd(1:3,3);
theta1 = atan2(p(2),p(1));
theta1(2,1) = theta1(1,1) - pi;

% ANGLES 2 & 3
% Find the rotation matrix associated with theta1, both sol'ns.
p1 = [q(:,4); 1];                   % At the intersection of axes 4,5,and 6
r2 = q(:,2);                        % Some point on the 2nd axis
r3 = q(:,3);                        % Some point on the 3rd axis
for i = 1:size(theta1);
    e1(:,:,i) = screw(w(:,1),v(:,1),theta1(i));
    g1(:,i) = inv(e1(:,:,i))*gd*inv(gst)*p1;
    % Here we've defined g1 = (e^-1)gd(gst^-1)p1, so that g1 = e2e3e4e5e6p1,
    %  but we've chosen p1 on axes 4,5,and 6; so e2e3p1 = g1
    
    % Now apply subproblem 2':
    [theta2, theta3] = sp2prime(w(:,2),w(:,3),p1(1:3),g1(1:3,i),r2,r3);
    
    % Angles 4 & 5
    p2 = [410 0 0 1]';              % A point on axis 6, not on 4 or 5
    r45 = [410 0 320]';             % The intersection of axes 4 and 5
    % Find all four versions of e2 and e3:
    for j = 1:size(theta2)
        e2(:,:,j) = screw(w(:,2),v(:,2),theta2(j));
        e3(:,:,j) = screw(w(:,3),v(:,3),theta3(j));
        g2(:,j) = inv(e3(:,:,j))*inv(e2(:,:,j))*inv(e1(:,:,i))...
            *gd*inv(gst)*p2;
        % Here we've found g2 = (e^-3)(e^-2)(e^-1)gd(gst^-1)p2 = e4e5e6p2
        %  but we choose p2 on axis 6 (and not on 4 or 5), so g2 = e4e5p2
        % The mod business in e^-1 ensures that theta1 is appropriately
        %  indexed to theta2 and theta3 (i=[1,2,3,4], j=[1,2,1,2])
        [theta4 theta5] = sp2(w(:,4),w(:,5),p2(1:3)-r45,g2(1:3,j)-r45);
        
        % Angle 6
        p3 = [150 0 260 1]';        % A point not on axis 6, for use in sp1
        r6 = [410 0 0]';            % A point on axis 6, for use in sp1
        for k = 1:size(theta5)
            e4(:,:,k) = screw(w(:,4),v(:,4),theta4(k));
            e5(:,:,k) = screw(w(:,5),v(:,5),theta5(k));
            g3(:,k) = inv(e5(:,:,k))*inv(e4(:,:,k))*inv(e3(:,:,j))...
                *inv(e2(:,:,j))*inv(e1(:,:,i))*gd*inv(gst)*p3;
            % Pick g3 = (e^-5)(e^-4)(e^-3)(e^-2)(e^-1)gd(gst^-1)p3 = e6p3
            %  so that g3 = e6p3.  Apply subproblem 1.
            theta6(k) = sp1(w(:,6),p3(1:3)-r6,g3(1:3,k)-r6);
            angles(i*j*k,:) = [theta1(i) theta2(j) theta3(j) ...
                theta4(k) theta5(k) theta6(k)];
            
            % Now test for joint limits and keep only the valid solutions:
            num = 4*mod(i+1,2) + 2*mod(j+1,2) + k;
            for m = 1:6
                if (angles(i*j*k,m) < jlimhigh(m) ...
                        & angles(i*j*k,m) > jlimlow(m))
                    FLAG = 1;
                else
                    fprintf('Beyond joint limits in sol''n ');
                    fprintf('%i: theta%i = %d\n',num, m, angles(i*j*k,m));
                    FLAG = 0;
                    break;
                end
            end
            if FLAG == 1
                theta(index,:) = angles(i*j*k,:);
                index = index + 1;
            end
            
        end
    end
end

offsets = [0 -pi/2 pi 0 -pi/2 0];
offsets = repmat(offsets,index-1,1);
theta
thetanew = theta + offsets