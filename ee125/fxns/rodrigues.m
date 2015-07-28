function R = rodrigues(omega, theta)
% RODRIGUES calculates the rotation matrix exp(omegahat*theta), given the
% axis of rotation omega and the angle of rotation theta, using Rodrigues's
% formula:
%       
%       R = I + w*sin(th) + w^2*(1-cos(th))
%
% where w is the skew symmetric matrix constructed from omega such that 
% cross(omega,v) = w*v (see SSM).  (Ordinarily, the matrix exponential is
% an infinite series, but if the matrix is skew-symmetric, as omegahat is,
% then the series can be written as the weighted sine and cosine fxns
% above.)

%-------------------------------------------------------------------------%
% Revised: 07/02/14
%   -vectorized input, tensorized output.
% Revised: 06/09/10
%   -edited description
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

% takes the matrix M, "turns it" into rows x slices, then stacks length(v)
% scaled versions of M along the columns of the resulting tensor, which has
% dimsions size(M,1) x length(v) x size(M,2)
matrixScale = @(M,v)(reshape(kron(v(:),M),[size(M,1),length(v),size(M,2)]));
omegahat = ssm(omega);

R = matrixScale(eye(length(omega)),ones(size(theta))) +...
    matrixScale(omegahat,sin(theta)) +...
    matrixScale(omegahat^2,1-cos(theta));

end