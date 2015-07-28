function EX = fa(D,p)
% FA    Factor analysis
%   FA takes data matrix D with *columns* as instances and a dimension p of
%   the latent variables and returns the expected value of the hidden units
%   given the data.
%
% note that this fxn assumes that each *column* of D is an instance
%
% We assume:
%   X ~ N(0,I)
%   Y ~ N(mu + lambda*x,Psi)
% Dimensions:
%           %-------------------------------------------%
%           %   x:          1Xp     |   M:          pXp %     
%           %   y:          1Xn     |   D:          nXN %
%           %   Lambda:     nXp     |   Y:          nXN %
%           %   Psi:        nXn     |   EX:         pXN %
%           %                       |   EXX_sum:    pXp %
%           %-------------------------------------------%
%-------------------------------------------------------------------------%
% Revised: 11/03/10
%   -corrected some mistakes!
% Cribbed: 11/03/10
%   -from your old cs281a assignment (see below)
% Adapted: 11/02/10
%   -from your old cs281_5.m
% Created: 11/16/04
%   by JGM
%-------------------------------------------------------------------------%

% init params
numiter = 100;
N = size(D,2);

% precondition
mu = sum(D,2)/N;                        % compute means
Y = D - repmat(mu,1,N);                 % de-mean

% init
L = rand(size(D,1),p);
P = rand(size(D,1));

for n = 1:numiter
    
    % E step
    Minv = eye(p) + L'/P*L;
    EX = Minv\L'/P*Y;
    EXX_sum = N*inv(Minv) + EX*EX';     % only need the *sum* of <xn*xn'>
%     M = eye(p) - L'/(P + L*L')*L;
%     EX = M*L'/P*Y;
%     EXX_sum = N*M + EX*EX';             % only need the *sum* of <xn*xn'>
    
    % M step
    Lold = L;
    L = Y*EX'/EXX_sum;
    if norm(Lold - L) < 0.00001;
        break
    end
    
    P = diag(diag(Y*Y' - L*EX*Y'))/N;   % the first diag removes the 
                                        %  diagonal components from its 
                                        %  matrix argument; the second 
                                        %  reconstructs a diagonal matrix 
                                        %  from these elements.
end

fprintf('Exited after %i iterations with error %d\n',n,norm(Lold - L));

end