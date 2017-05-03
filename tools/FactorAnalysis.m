function [M,b,SigmaYX,W,SigmaXY] = FactorAnalysis(Y,Nltnt)
% FactorAnalysis    Factor analysis
%
%   FactorAnalysis implements an EM algorithm to learn a linear-Gaussian
%   generative model for data Y:
%
%       p(x)    = N(0,I),
%       p(y|x)  = N(M*x + b, SigmaYX),      SigmaYX diagonal!
%   =>  p(x|y)  = N(W*(y-b), SigmaXY)
%
%   That is, the function learns, and returns, the generative parameters M,
%   b, and SigmaYX.  It also returns the parameters of the posterior, W and
%   SigmaXY.  The dimensionality of X, Nltnt, must also be provided.
%
%   N.B.: This funtion expects the *columns* of Y to be samples!

%-------------------------------------------------------------------------%
% Revised: 04/12/16
%   -renamed from fa.m to FactorAnalysis.m
%   -largely rewrote:
%       --made more efficient; 
%       --renamed all variables;
%       --implemented cross-entropy calculation, and hooked exit to it
%       --changed output variables
%       --functionized
% Revised: 11/03/10
%   -corrected some mistakes!
% Cribbed: 11/03/10
%   -from your old cs281a assignment (see below)
% Adapted: 11/02/10
%   -from your old cs281_5.m
% Created: 11/16/04
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nitermax = 1000;
[Nobsv,Nsamples] = size(Y);
thr = 0.00001;     % tiny!


% de-mean the observations
b = mean(Y,2);
Y = Y - b;
YY = Y*Y'/Nsamples;

% init params
M = rand(Nobsv,Nltnt);
SigmaYX = rand(Nobsv,Nobsv);
XNtrpY = Inf;
CONT = 1;
iIter = 1;

% loop over E step, M step, cross-entropy calculation
while (iIter < Nitermax)&&CONT
    [XX,YX,W,SigmaXY] = getXpctSuffStats(YY,M,SigmaYX);
    [M,SigmaYX] = getMLparams(XX,YX,YY);
    [XNtrpY,CONT] = getCrossEntropy(Y,XNtrpY,M,SigmaYX,W,thr);
    iIter = iIter + 1;
end
fprintf('Exited after %i iterations\n',iIter);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [XX,YX,W,SigmaXY] = getXpctSuffStats(YY,M,SigmaYX)
% E step: compute expected suff. stats. of posterior, given these params

%%% this is equivalent, but let's not take so many inverses (see labnotes)
% InfoXY = eye(size(M,2)) + M'/SigmaYX*M;
% W = InfoXY\M'/SigmaYX;              % this is also the "Kalman gain"
% XX = inv(InfoXY) + W*YY*W';
W = M'/(M*M' + SigmaYX);
SigmaXY = eye(size(W,1)) - W*M;
XX = SigmaXY + W*YY*W';
YX = YY*W';


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [M,SigmaYX] = getMLparams(XX,YX,YY)
% M step: solve for the best params, given these suff. stats.

M = YX/XX;
SigmaYX = YY - M*YX';
SigmaYX = diag(diag(chol(SigmaYX*SigmaYX')));      % for numerical stability

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [XNtrpYnew,CONT] = getCrossEntropy(Y,XNtrpYold,M,SigmaYX,W,thr)
% compute cross entropy, and issue break if it stops decreasing

XNtrpYnew = GaussianCrossEntropy(Y,M*M' + SigmaYX,'cvrn',1);
fractionalImprovement = (XNtrpYold - XNtrpYnew)/abs(XNtrpYold);
if (fractionalImprovement < thr)
    fprintf('exiting with delta cross-entropy below tolerance\n');
    CONT = 0;
else
    CONT = 1;
end

end
%-------------------------------------------------------------------------%