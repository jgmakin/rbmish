function [beta,Rsq,Yres,Yhat,p] = linregress(X,Y,varargin)
% linregress    Linear regression
%
% USAGE:
%   beta = linregress(X,Y)
%   [beta,Rsq,Yres,Yhat,p] = linregress(X,Y)
%   [beta,RsqCV,YresCV,YhatCV,pCV] = linregress(X,Y,'LOO')
%   [beta,RsqCV,YresCV,YhatCV,pCV] = linregress(X,Y,'LOO',lambda)
% 
% linregress(X,Y) regresses the columns of Y onto X; that is, it finds the
% most likely coefficients beta under the model:
%
%   Y ~ N(X*beta,Sigma)
%
% Both X and Y can be matrices, whose *rows* are "samples" or "trials."  
% Hence:
%
%   size(X)     = Nsamples x Nfeatues
%   size(Y)     = Nsamples x Noutputvars
%   size(beta)  = Nfeatures x Noutputvars
% 
% One can also request the corresponding coefficients of determination,
% Rsq; residuals, Yres; model estimates, Yhat; and p-values (under the
% f-test) for fit quality.
%
% If 'LOO' is provided as the third argument, the reported Rsq, Yres, Yhat 
% and p values are leave-out-out cross-validated.  But NB that the returned
% coefficients beta are *not* cross-validated (this would require returning
% Nsamples different (Nfeatures-1) x Noutputvars matrices.
%
% The optional fourth argument is the L2-regularization penalty (otherwise 
% assumed 0).
% 
% NB that, unlike matlab's version, the column of ones is NOT included by
% default.
% 
% For a derivation of the loopless computation of residuals under LOOCV,
% see your lab notes.

%-------------------------------------------------------------------------%
% Version history for linrgsLOOCV.m:
%
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
% Revised: 08/17/16
%   -merged with linrgsLOOCV
% Revised: 06/06/16
%   -added p-value calculation
% Created: 03/22/12
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nsamples = size(Y,1);

% preliminaries
crossvalidation = 'none';
if ~isempty(varargin), crossvalidation = varargin{1}; end

% L2 regularization?
if length(varargin)>1
    lambda = varargin{2};
    Xpinv = (X'*X + lambda*eye(size(X,2)))\X';
else
    Xpinv = (X'*X)\X';
    %%%% replace with computation based on qr decomposition
end

% parameters (MLE under assumed normal distribution p(y|x;beta))
beta = Xpinv*Y;



% all the other outputs
if nargout > 1
    
    % residuals
    switch crossvalidation
        case 'LOO'
            % cross-validated residuals
            if Nsamples < 1000
                H = X*Xpinv;                    % hat matrix
                M = eye(Nsamples) - H;          % residual matrix
                Yres = (M*Y)./diag(M);
                %%% Yres = diag(1./diag(B))*B*Y; % *much* slower
            else
                ResAll = Y - X*beta;
                Yres = zeros(size(ResAll),'like',X);
                for i = 1:Nsamples
                    Yres(i,:) = ResAll(i,:)/(1 - X(i,:)*Xpinv(:,i));
                    % cross-validated residuals
                end
                %%% this is faster than arrayfun
            end
            % estimate of Y
            if nargout > 3, Yhat = Y - Yres; end
        otherwise
            Yhat = X*beta;
            Yres = Y - Yhat;
    end
    
    % coefficient of determination
    SStot = sum((Y - mean(Y,1)).^2);
    SSerr = sum(Yres.^2,1);
    Rsq = 1 - SSerr./SStot;
    
    % p value in comparison to a "model" that uses only the mean
    if nargout > 4
        Nvars = size(X,2);
        f = ((Nsamples-Nvars)/(Nvars-1))*((SStot-SSerr)./SSerr);
        p = 1 - fcdf(f,Nvars-1,Nsamples-Nvars);
    end
    
end

end
