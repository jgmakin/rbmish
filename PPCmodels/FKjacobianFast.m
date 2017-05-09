function J = FKjacobianFast(th,roboparams)
% FKjacobianFast    Fast computation of the forward-kinematics jacobian
%
%   USAGE:
%       J = FKjacobianFast(th,roboparams)
%
%   FKjacobianFast takes a matrix of joint angles, th, of size (Nexamples x
%   Ndims), and the parameters structure; and returns a tensor (Nexamples x
%   Ndims x Ndims) of Jacobians of the forward kinematics, J, at those
%   points.
%
%   NB the dimensions!!

%-------------------------------------------------------------------------%
% Revised: 07/07/14
%   -transposed the expected input (th) and output
% Revised: 07/04/14
%   -matricized input pos, tensorized output J
% Adapted: 06/03/11
%   -from marginalErrors.m
%-------------------------------------------------------------------------%

Ndims = length(roboparams.armlengths);

if Ndims == 1
    J = -roboparams.armlengths(1)*sin(th);
else
    % init params
    L1 = roboparams.armlengths(1);
    L2 = roboparams.armlengths(2);
    
    thSum = sum(th,2);
    J(:,1,1) = -L1*cos(th(:,1)) - L2*cos(thSum);
    J(:,1,2) = -L2*cos(thSum);
    J(:,2,1) = -L1*sin(th(:,1)) - L2*sin(thSum);
    J(:,2,2) = -L2*sin(thSum);
end

end