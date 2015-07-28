function CvrnMX = covX2covMX(CvrnX,M)
% A strange function, perhaps. CvrnXhat is a tensor of covariance matrices:
%   
%   Ncols x Ncols x T
%
% where Ncols = Nstates or Noutputs (or etc.?).  Naively, mutiplying left 
% and right sides of each covariance matrix by M requires a loop through T 
% (generally long).  Instead, we use two nested loops of length Ncols
% (generally small).
%
% For size(CvrnX) = 4x4, M = 2x4, T = 1000, this function is indeed orders
% of magnitude faster than its naive equivalent.
%
% NB: This also turns out to be much faster than the tensor operation!!
%
%   Mtensor = repmat(M,[1,1,T]);
%   MSigma = tensorOp(permute(Mtensor,[1,3,2]),permute(CvrnX,[1,3,2]));
%   CvrnMX = permute(tensorOp(MSigma,permute(Mtensor,[2,3,1])),[1,3,2]);
%
% presumably b/c of the reshape operations (?!).  So don't be tempted to
% vectorize this.

%-------------------------------------------------------------------------%
% Revised: 07/15/14
%   -compared efficiency with tensorOp
% Created: ??/??/14
%   bu JGM
%-------------------------------------------------------------------------%



% Ns
T = size(CvrnX,3);
[Nrows,Ncols] = size(M);

% malloc
foo = zeros(Nrows^2,Ncols^2,'like',CvrnX);

% loop
for i = 1:Nrows                               % get vectorized outer
    for j = 1:Nrows                           %   products of all the 
        ind = (i-1)*Nrows + j;                %   rows of M
        op = M(i,:)'*M(j,:);
        foo(ind,:) = op(:);
    end
end
caca = foo*reshape(CvrnX,[Ncols^2,T]);
CvrnMX = reshape(caca,[Nrows,Nrows,T]);

end

