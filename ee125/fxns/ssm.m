function omegahat = ssm(omega)
% SSM Skew Symmetric Matrix
%   Creates a skew symmetric matrix omegahat of the form
%   
%   [0 -w3 w2; w3 0 -w1; -w2 w1 0]
%   
%   given a vector omega, such that omega x v = omegahat*v. 
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%

omegahat = [0 -omega(3) omega(2);...
            omega(3) 0 -omega(1);...
           -omega(2) omega(1) 0];