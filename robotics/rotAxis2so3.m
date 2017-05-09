function omegahat = rotAxis2so3(omega)
% rotAxis2so3  Convert rotation axis to equivalent skew symmetric matrix
%
% omegahat = rotAxis2so3(w) converts the rotation axis w into the matrix
% omegahat \in so(3):
%   
%   [0 -w3 w2; w3 0 -w1; -w2 w1 0],
%
% so that omega x v = omegahat*v.
%
% NB: omega must have size 3 x Nexamples, and the output omegahat has size
% 3 x 3 x Nexamples.

%-------------------------------------------------------------------------%
% Revised: 11/20/16
%   -vectorized to work on more than one rotation axis at once 
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

% a vector of zeros
zz = zeros(1,1,size(omega,2));

% each slice is a skew-symmetic matrix, for Nexamples slices
omega = shiftdim(omega,-1);
omegahat = [zz -omega(1,3,:) omega(1,2,:);...
            omega(1,3,:) zz -omega(1,1,:);...
           -omega(1,2,:) omega(1,1,:) zz];
       

end