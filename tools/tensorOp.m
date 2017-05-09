function NQ = tensorOp(NP,PQ)
% tensorOp      Multiply many matrices by many matrices
%
% USAGE:
%   NQ = tensorOp(NP,PQ)
%
% Given:
%
%   (n x p x m) tensor NP
%
% *where m is the number of matrices*, and
%
%   (p x q x m) tensor PQ
%
% compute all m of the matrix multiplications, and store the result in
%
%   (n x q x m) tensor NQ.
%
% Notice that to multiply the set of matrices in NP by a set of *vectors*
% in PQ, PQ must have size (p x 1 x m).
%
% To multiply a set of matrices NP by a *single vector* P, use:
%
%   N = tensorOp(NP,repmat(P,[1,1,m])

%-------------------------------------------------------------------------%
% Revised: 11/18/16
%   -changed the ordering of inputs and outputs!!  Now the number of
%   matrices, m, is the *final* (third), rather than second, dimension.
% Created: 07/04/14
%   by JGM
%-------------------------------------------------------------------------%

NQ = sum(permute(NP,[1,4,3,2]).*permute(PQ,[4,2,3,1]),4);

end