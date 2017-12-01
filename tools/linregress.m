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
%
% Variable arguments:
%
%   'cross validate':   'LOO' | 'none'
%
% If 'LOO', the reported Rsq, Yres, Yhat and p values are leave-out-out 
% cross-validated.  But NB that the returned *coefficients* beta are *not*
% cross-validated (this would require returning Nsamples different 
% (Nfeatures-1) x Noutputvars matrices.
%
%   'L2 regularize'     <value>
%
% If <value> is greater than 0, it will be interpreted as lambda, the L2
% regularizer penalty.
%
%   'pad with ones'     <value>
%
% NB that, unlike matlab's version, the column of ones is NOT included by
% default.  But if <value> is not 0, it will.  NB!!  The column of ones is
% always the *last* column.
% 
% For a derivation of the loopless computation of residuals under LOOCV,
% see your lab notes.

%-------------------------------------------------------------------------%
% Version history for linrgsLOOCV.m:
% Revised: 10/10/17
%   -changed variable arguments to use your defaulter function.
%   -added 'pad with ones' as a variable input argument.
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



% variable arguments
crossvalidation = defaulter('cross validate','none',varargin{:});
lambda_L2 = defaulter('L2 regularize',0,varargin{:});
pad_with_ones = defaulter('pad with ones',0,varargin{:});


% Ns
Nsamples = size(Y,1);

% note that we always put the ones in the *last* column
if pad_with_ones, X = [X, ones(size(X,1),1)]; end

% L2 regularization?
if lambda_L2 > 0
    Xpinv = (X'*X + lambda_L2*eye(size(X,2)))\X';
    %%%% get rid of construction of identity matrix....
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
