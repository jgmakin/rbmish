function X = dareJGM(F,G1,G2,H)
% USAGE:
%
%   P = dareJGM(A,B,R,C'*Q*C)
%
% (For some inconceivable reason, Laub uses the other set of letters.)
%
% If you don't have Matlab's control toolbox---or you've used up all the
% licenses---you will need to use this function.
%
% This implementation follows Laub 1979, "A Schur Method for Solving 
% Algebraic Riccati Equations."

%-------------------------------------------------------------------------%
% Created: 04/03/14
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nstates = size(F,1);

% define
G = G1*inv(G2)*G1';
M = inv(F)';
Z = [F + G*M*H, -G*M; -M*H, M];

% get the Schur decomposition; that is, get the quasi-upper-triangular 
% Schur matrix T and the unitary matrix U s.t. U*T*U' = Z, and (therefore) 
% U'*Z*U = T.  
% Using 'real', we also ensure that the 2x2 blocks on the block diagonal of
% T correspond to the complex eigenvalues, while the 1x1 blocks correspond
% to the real eigenvalues (i.e., Murnaghan-Wintner canonical form).
[Ut,T] = schur(Z,'real');

% now reorder the blocks so that the stable eigenvalues (within the unit
% circle) are all in the upper-left block, and the unstable all in the
% lower-right block.  Here, U*S*U' = Z still, and thus U'*Z*U = S.
[U,S] = ordschur(Ut,T,'udi');

% grab the upper-left and the lower-left blocks of U 
U_11 = U(1:Nstates,1:Nstates);
U_21 = U(Nstates+1:end,1:Nstates);

% the solution
X = U_21/U_11;

%%% check
% F'*X*F - X - F'*X*G1*inv(G2 + G1'*X*G1)*G1'*X*F + H


end







