
clear;clc;
gd = [-0.4317 -0.5344 0.7266 325.3; 0.5152 -0.8073 -0.2877 92.3;...
    0.7404 0.2502 0.6239 508.6; 0 0 0 1];
gst = [eye(3) [410; 0; 410]; 0 0 0 1];
theta1 = 0.2765;
e1 = screw([0 0 1]',[0 0 0]',theta1);
p1 = [410 0 320 1]';
g1 = inv(e1)*gd*inv(gst)*p1;
q = g1(1:3);
p = p1(1:3);
r = [150 0 260]';
s = [150 0 0]';
u2 = p - r;
u1 = q - s;
norm(u2);
norm(u1);

c = [361 0 423]';
for i=1:2
    v1 = c - s;
    v2 = c - r;
    theta2 = sp1([0 1 0]',u1,v1);
    theta3 = sp1([0 1 0]',u2,v2);
    fprintf('angle 2 = %d\n', theta2);
    fprintf('angle 3 = %d\n', theta3);
    e2 = screw([0 1 0]',cross(-[0 1 0]',s),theta2);
    e3 = screw([0 1 0]',cross(-[0 1 0]',r),theta3);
    p2 = [410 0 410 1]';
    g2 = inv(e3)*inv(e2)*inv(e1)*gd*inv(gst)*p2;
    [theta41 theta51 theta42 theta52] =....
        sp2([1 0 0]',[0 1 0]',p2(1:3),g2(1:3),p1(1:3))
    
    
    c = [-61 0 423]';
end