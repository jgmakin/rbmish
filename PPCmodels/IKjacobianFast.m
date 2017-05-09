function J = IKjacobianFast(pos,roboparams)
% IKjacobianFast    Fast Jacobian of the inverse kinematics
%
%   USAGE:
%       IKjacobianFast(pos,roboparams)
% 
%   IKjacobianFast takes a matrix of end-effector positions, pos, of size
%   (Nexamples x Ndims), and the parameters structure; and returns a tensor
%   (Nexamples x Ndims x Ndims) of Jacobians of the inverse kinematics, J, 
%   at those points.
%
%   NB the dimensions!!

%-------------------------------------------------------------------------%
% Revised: 07/07/14
%   -transposed input (pos) and output
% Revised: 07/04/14
%   -matricized input pos, tensorized output J
% Adapted: 05/31/11
%   -from marginalErrors.m
%   -This has been checked in a dozen different ways, so the math is surely
%   right.
%   by JGM
%-------------------------------------------------------------------------%


Ndims = length(roboparams.armlengths);

if Ndims == 1
    J = -(roboparams.armlengths(1)^2 - pos.^2).^(-1/2);
else
    
    L1 = roboparams.armlengths(1);
    L2 = roboparams.armlengths(2);
    
    rsq = pos(:,1).^2 + pos(:,2).^2;
    gamma = (L1^2 + L2^2 - rsq)/(2*L1*L2);
    delta = (rsq + L1^2 - L2^2)./(2*L1*sqrt(rsq));
    
    J(:,1,1)= -pos(:,2)./rsq + 1./sqrt(1-delta.^2).*pos(:,1).*...
        (sqrt(rsq)-L1*delta)/L1./rsq;
    J(:,1,2) = pos(:,1)./rsq + 1./sqrt(1-delta.^2).*pos(:,2).*...
        (sqrt(rsq)-L1*delta)/L1./rsq;
    J(:,2,1) = -pos(:,1)./(sqrt(1-gamma.^2)*L1*L2);
    J(:,2,2) = -pos(:,2)./(sqrt(1-gamma.^2)*L1*L2);
end

end