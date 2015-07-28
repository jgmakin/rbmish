function [beta,Yhat,Rsq] = linregress(X,Y)
% Yes, you wrote your own, for no obvious reason

%-------------------------------------------------------------------------%
% Created: 03/22/12
%   by JGM
%-------------------------------------------------------------------------%


beta = (X'*X)\X'*Y;
Yhat = X*beta;

Ybar = mean(Y,1);
Res = Y - Yhat;
SSerr = sum(Res.^2,1);
SStot = sum((Y - repmat(Ybar,size(Y,1),1)).^2);
Rsq = 1 - SSerr./SStot;


end