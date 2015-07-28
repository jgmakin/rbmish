function [theta1, theta2] = sp2(w1,w2,u,v);
% SP2 computes the inverse kinematics Paden-Kahan subproblem 2.
%  Given two axes of rotation w1 and w2, an initial point p and the final
%  point q, and the intersection point of the two axes, r, (see Fig. 3.8 in
%  Murray, Lee, Sastry), it returns the values of theta1 and theta2, the
%  angles of rotation for the two axes.

alpha = ((w1'*w2)*w2'*u - w1'*v)/((w1'*w2)^2-1);
beta  = ((w1'*w2)*w1'*v - w2'*u)/((w1'*w2)^2-1);
gammasquared = (norm(u)^2 - alpha^2 - beta^2 - 2*alpha*beta*w1'*w2)/...
    (norm(cross(w1,w2))^2);
if gammasquared < 0                         % These crazy results will be 
    theta1 = 20000;                         %  produced just in the case 
    theta1(2,1) = 20000;                    %  that there are no sol'ns to
    theta2 = 20000;                         %  the subproblem.  They'll be
    theta2(2,1) = 20000;                    %  caught and discarded by the 
else                                        %  joint-limits check @ the end
    gamma = sqrt(gammasquared);
     
    z = alpha*w1 + beta*w2 + gamma*(cross(w1,w2));
    theta1 = sp1(-w1,v,z);
    theta2 = sp1(w2,u,z);
    
    z = alpha*w1 + beta*w2 - gamma*(cross(w1,w2));
    theta1(2,1) = sp1(-w1,v,z);
    theta2(2,1)= sp1(w2,u,z);
end