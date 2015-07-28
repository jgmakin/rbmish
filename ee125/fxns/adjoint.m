function Adg = adjoint(g)
% ADJOINT   Adjoint tranformation for twist transformations.
%   ADJOINT(G) returns the 6x6 adjoint transformation that transforms
%   twists (members of SE(3)) from one coordinate frame to another.

R = g(1:3,1:3);
p = g(1:3,4);
phat = ssm(p);
Adg = [R phat*R; zeros(3) R];