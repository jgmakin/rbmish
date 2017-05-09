function R = rodrigues(omega,th)
% Rodrigues's rotation formula
%
% USAGE:
%   R = rodrigues(omega,th)
% 
% Calculates the rotation matrix exp(omegahat*th), given the axis of
% rotation omega and the angle of rotation th, using Rodrigues's formula:
%       
%       R = I + omegahat*sin(th) + omegahat^2*(1-cos(th))
%
% where omegahat is the skew symmetric matrix constructed from omega such 
% that cross(omega,v) = w*v (see SSM).  (Ordinarily, the matrix exponential
% is an infinite series, but if the matrix is skew-symmetric, as omegahat 
% is, then the series can be written as the weighted sine and cosine fxns
% above.)
%
% Th angle of rotation th is a vector of Nexamples elements, but the axis 
% of rotation omega can have size *either* (3 x 1) or (3 x Nexamples).

%-------------------------------------------------------------------------%
% Revised: 11/20/16
%   -rewrote to take advantage of the new rotAxis2so3.m (formerly ssm.m),
%   which is now tensorized (can operate on multiple axes at once).  This
%   function now works for many axes omega (3 x Nexamples)--as long as the
%   number of columns matches the the number of elements in th.
% Revised: 11/18/16
%   -changed output ordering from (3 x Nexamples x 3) to (3x3 x Nexamples),
%   in accordance with changes to tensorOp.m.
% Revised: 07/02/14
%   -vectorized input, tensorized output.
% Revised: 06/09/10
%   -edited description
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

%%%matrixScale = @(M,v)(reshape(kron(v(:),M),[size(M,1),length(v),size(M,2)]));
Nexamples = length(th);
matrixScale = @(M,v)(M.*shiftdim(v(:),-2));
omegahat = rotAxis2so3(omega);
omegahatSq = tensorOp(omegahat,omegahat);

R = repmat(eye(3),[1,1,Nexamples]) +...
    matrixScale(omegahat,sin(th)) +...
    matrixScale(omegahatSq,1-cos(th));

end
































