function params = refitTransitionNoise(X,endinds,params)
% "fix" the parameters from the fit of the dynamical mode

%-------------------------------------------------------------------------%
% Revised: 03/11/15
%   -added new case, which is now the default, to refit the whole A matrix
%   on the subsampled trajectories, rather than just taking A^m.  Some spot
%   checks show this to be superior (in terms of LOOCVed R^2).
% Created: ??/??/13
%   by JGM
%-------------------------------------------------------------------------%

% params
m = params.m;


% prepare regressions
XFuture = X(~ismember(1:size(X,1),endinds+1),:);
XPast = X(~ismember(1:size(X,1),endinds),:);
Xpinv = (XPast'*XPast)\XPast';

switch params.BINMETHOD
    case 'slidingwindow'
        A = params.A;
    case 'oldwayofdoingit'
        A = params.A^m;
    otherwise
        fprintf('Refitting all of A on subsampled trajectories ')
        fprintf('(you discovered this to be the best)\n\n');
        A = (Xpinv*XFuture)';
end


% regression residuals with LOOCV
Res = XFuture - XPast*A';
XResCV = zeros(size(Res));
for i = 1:size(XPast,1)
    XResCV(i,:) = Res(i,:)/(1 - XPast(i,:)*Xpinv(:,i));
end
muX = mean(XResCV);
SigmaX = cov(XResCV);


% how well did you do?
SSerr = sum(XResCV.^2);
SStot = sum((XFuture - repmat(mean(XFuture,1),size(XFuture,1),1)).^2);
RsqCV = 1 - SSerr./SStot;
fprintf('A fit with R^2 = %0.3f\n\n',RsqCV);


% store
params.A = A;
params.SigmaX = SigmaX;
params.muX = muX';


end