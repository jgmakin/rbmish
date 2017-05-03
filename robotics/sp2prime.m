function [theta1, theta2] = sp2prime(w1,w2,p,q,r1,r2) 
%SP2PRIME computes the angles of rotation theta1 and theta2
% for two parallel axes, using Paden-Kahan subproblem 2'. Given two axes of
% rotation, w1 and w2, an initial point p and a final point q, and two
% points on the respective axes, r1 and r2, it returns the angles of
% rotation for each axis.
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

ce1 = r1 + w1*w1'*(q-r1);
ce2 = r2 + w2*w2'*(p-r2);
R1 = norm(cross(q-r1,w1));
R2 = norm(cross(p-r2,w2));
D = norm(ce2-ce1);
d = (R1^2-R2^2+D^2)/(2*D);
s = (ce2-ce1)/D;
t = cross(s,w1);
if R1^2 < d^2
    theta1 = 10000;                             % These spurious results 
    theta1(2,1) = 10000;                        %  will be caught by the 
    theta2 = 10000;                             %  joint limits check at 
    theta2(2,1) = 10000;                        %  the end
else
    c = ce1 + d*s + sqrt(R1^2-d^2)*t;           % First set of sol'ns
    theta1 = sp1(-w1,q-r1,c-r1);
    theta2 = sp1(w2,p-r2,c-r2);     
    
    c = ce1 + d*s - sqrt(R1^2-d^2)*t;           % Second set of sol'ns                                    
    theta1(2,1) = sp1(-w1,q-r1,c-r1);
    theta2(2,1) = sp1(w2,p-r2,c-r2);   
end