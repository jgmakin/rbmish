function varargout = linrgsLOOCV(X,Y,varargin)
% leave-one-out cross validation (LOOCV).
% USAGE:
%       [beta, RsqCV] = linrgsLOOCV(Xavg,Z);
%       [beta, RsqCV] = linrgsLOOCV(Xavg,Z,lambda);
%       [beta, RsqCV, YRes] = linrgsLOOCV(Xavg,Z);
%
% For a derivation, see your lab notes
%
% NB: does not include the column of ones!!  You have to add that.

%-------------------------------------------------------------------------%
% Revised: ??/??/??
%   -added if statement to use different methods depending on the number of
%   data vectors (rows of Y)
% Revised: 05/15/13
%   -added varargin for an L2 penalty
%   -replaced the calculation of the "hat matrix" with one that requires
%   far fewer operations---but uses a loop
% Revised: 05/??/13
%   -added varargout for the residuals
% Created: 04/01/12
%   by JGM
%-------------------------------------------------------------------------%

N = size(Y,1);

% compute the right pseudoinverse
if isempty(varargin)
    Xpinv = (X'*X)\X';
else
    lambda = varargin{1};                               % L2 penalty
    Xpinv = (X'*X + lambda*eye(size(X,2)))\X';
end
beta = Xpinv*Y;

% compute the residuals
if N < 1000
    H = X*Xpinv;                    % hat matrix
    B = (eye(N) - H);               % residual matrix
    ResCV = diag(1./diag(B))*B*Y;   % residuals
else
    Res = Y - X*beta;
    ResCV = zeros(size(Res),'like',X);
    for i = 1:size(X,1)
        ResCV(i,:) = Res(i,:)/(1 - X(i,:)*Xpinv(:,i));
    end
end

% compute the the cross-validated coefficient of determination
SSerr = sum(ResCV.^2);
Ybar = mean(Y,1);
SStot = sum((Y - repmat(Ybar,N,1)).^2);
RsqCV = 1 - SSerr./SStot;

% [beta, RsqCV]
varargout{1} = beta;
varargout{2} = RsqCV;
varargout{3} = ResCV;

end



    
    
    
    
    
    
    
    
    
    
    
    