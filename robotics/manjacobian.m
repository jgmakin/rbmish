function J = manjacobian(w,q,th)
% MANJACOBIAN   Manipulator spatial jacobian
%   This program calculates the spatial Jacobian for the AdeptSix 300 robot 
%   arm, given a joint configuration
%-------------------------------------------------------------------------%
% Revised: 10/08/10
%   -preallocated matrices and vectors "for speed"
% Adapted: 08/19/10
%   -from manjacobian in lab3 files (see for details)
%-------------------------------------------------------------------------%

% check
m = length(th);
if size(w,2) ~= m
    error('mismatch b/n number of axes and angles - jgm')
end
if size(q,2) ~= m
    error('mismatch b/n number of axis points and angles - jgm')
end
     
% malloc
v = zeros(3,m);
xi = zeros(6,m);
SE3 = zeros(4,4,m);
J = zeros(6,m);

% do 
for i = 1:m
    v(:,i) = cross(-w(:,i),q(:,i));
    xi(:,i) = [v(:,i); w(:,i)];
    SE3(:,:,i) = screw(w(:,i),v(:,i),th(i));
    g = eye(4);
    for j = 1:i-1;
        g = g*SE3(:,:,j);
    end
    J(:,i) = adjoint(g)*xi(:,i);
end

end