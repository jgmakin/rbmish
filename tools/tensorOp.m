function NQ = tensorOp(NP,PQ)
% tensorOp      Multiply many matrices by many matrices
%
% USAGE:
%   NQ = tensorOp(NP,PQ)
%
% Given:
%
%   (n x m x p) tensor NP
%
% *where m is the number of matrices*, and
%
%   (p x m x q) tensor PQ
%
% compute all m of the matrix multiplications, and store the result in
%
%   (n x m x q) tensor NQ.
%
% Notice that this still works when PQ is a matrix of size (p x m), i.e.
% when multiplying the matrices in NP by the set of *vectors* in PQ.
%
% To multiply a set of matrices NP by a *single vector* P, use:
%
%   N = tensorOp(NP,repmat(P,[1,m])

%-------------------------------------------------------------------------%
% Created: 07/04/04
%   by JGM
%-------------------------------------------------------------------------%

[n,m,p] = size(NP);
q = size(PQ,3);
Qfour(1,1:m,1:q,1:p) = permute(PQ,[2 3 1]);     % 1 x m x q x p
Tfour(1:n,1:m,1,1:p) = NP;                      % n x m x 1 x p
    
NQ = sum(Tfour(:,:,ones(q,1),:).*Qfour(ones(n,1),:,:,:),4); % n x m x q

end
