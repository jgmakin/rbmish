function R = expCoords2SO3(omegath)
% expCoords2SO3     Convert from "exponential coordianates" to SO(3)
%
% USAGE:
%   R = expCoords2SO3(omegath);
%
% We sometimes (e.g., the mocap data set) have "exponential coordinates" of
% a rotation, i.e. omega*th, with omega the axis, and th the amount, of the
% rotation.  We'd like to convert this product to a rotation matrix.
%
% The input must have size 3 x Nexamples (but Nexamples may be 1) 
%
% See also the inverse function, SO32expCoords.m.

%-------------------------------------------------------------------------%
% Created: 11/18/16
%   by JGM
%-------------------------------------------------------------------------%

% special case: no rotation
iNoRot = sum(abs(omegath)) < 3*eps;

theta = sqrt(sum(omegath.*omegath));
omega = omegath./theta;
omega(:,iNoRot) = repmat([1;0;0],[1,sum(iNoRot)]);  % arbitrary
R = rodrigues(omega,theta);
 

end