function [xi,XpctYX,CvrnYX,qXgivenY] = GaussianMixtureModel(Y,Nclasses,varargin)
% GaussianMixtureModel  Fit a GMM with expectation-maximization
%
% USAGE:
%   [xi,XpctYX,CvrnYX] = GaussianMixtureModel(Y,Nclasses);
%
% GaussianMixtureModel implements an EM algorithm to learn the following 
%   generative model for data Y:
%
%       p(X)    = Cat[xi]
%       p(Y|X)  = N(mu(x), Sigma(x))
%
%   That is, the function learns, and returns, the generative parameters 
%   xi, XpctYX = cat(mu{:}), CvrnYX = cat(Sigma{:}).  It also returns the
%   posterior class assignments, qXgivenY.  In addition to the data Y, the 
%   number of classes Nclasses must also be provided as an input.
%
%   N.B.: This funtion expects the *columns* of Y to be samples!

%-------------------------------------------------------------------------%
% Created: 02/09/17
%   by JGM
%-------------------------------------------------------------------------%

%%%%% TO DO
% (1) change to use log probabilities (or add as case)


% Booleans
TOANIMATE  = defaulter('animate',0,varargin{:});
%%% USELOGS = defaulter('use logs',1,varargin{:});

% Ns
Nitermax = 1000;
Nobsv = size(Y,1);
thr = 0.00001;     % tiny!

% init model params with k means
XpctYX  = zeros(Nobsv,Nclasses,'like',Y);
CvrnYX = zeros(Nobsv,Nobsv,Nclasses,'like',Y);
IDs = kmeans(Y', Nclasses);
xi = mean(IDs == (1:Nclasses));
for iClass = 1:Nclasses
    theseY = Y(:,IDs == iClass);
    XpctYX(:,iClass)   = mean(theseY,2);
    CvrnYX(:,:,iClass) = cov(theseY');
end

% initialize non-model vars
qXgivenY = permute(IDs == (1:Nclasses),[3,1,2]);
XNtrpY = Inf;
CONT = 1;
iIter = 1;

% loop over E step, M step, cross-entropy calculation
while (iIter < Nitermax)&&CONT
    if TOANIMATE, plotGMMfit(Y,xi,XpctYX,CvrnYX,qXgivenY); end
    [XNtrpYnew,WtdXpctY,WtdXpctYY,qXgivenY,qpX] =...
        getXpctSuffStats(Y,xi,XpctYX,CvrnYX);
    [xi,XpctYX,CvrnYX] = getMLparams(qpX,WtdXpctY,WtdXpctYY);
    [XNtrpY,CONT] = getCrossEntropy(XNtrpY,XNtrpYnew,thr);
    fprintf('cross entropy: %.3f\n',XNtrpY); 
    iIter = iIter + 1;
end
fprintf('Exited after %i iterations\n',iIter);

qXgivenY = permute(qXgivenY,[3,2,1]);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [XNtrpY,WtdXpctY,WtdXpctYY,qXgivenY,qpX] =...
    getXpctSuffStats(Y,xi,XpctYX,CvrnYX)
% Most variables are in (Ndims x Nsamples x Nclasses) format

% Ns
Nsamples = size(Y,2);
Nclasses = length(xi);

% joint--i.e., *unnormalized*-wrt-X (sample-based) posterior dstrb
qXY = zeros(1,Nsamples,Nclasses,'like',Y);
for iClass = 1:Nclasses
    qXY(:,:,iClass) = xi(iClass)*mvnpdf(Y',XpctYX(:,iClass)',gather(CvrnYX(:,:,iClass)))';
end
qy          = sum(qXY,3);
qXgivenY    = qXY./qy;
YqXgivenY   = Y.*qXgivenY;

% expected sufficient statistics
XNtrpY      = -sum(log(qy))/Nsamples;
qpX         = sum(qXgivenY,2)/Nsamples;
WtdXpctY    = sum(YqXgivenY,2)/Nsamples;
WtdXpctYY   = sum(permute(Y,[1,3,2]).*permute(YqXgivenY,[4,1,2,3]),3)/Nsamples;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [xi,XpctYX,CvrnYX] = getMLparams(qpX,WtdXpctY,WtdXpctYY)

% M step: maximum-likelihood parameters
xi = permute(qpX,[2,3,1]);
mu = WtdXpctY./qpX;
XpctYX = permute(mu,[1,3,2]);
CvrnYX = permute(WtdXpctYY,[1,2,4,3])./qpX - mu.*permute(mu,[2,1,3]);

% unfortunate addendnum for numerical stability
for iClass = 1:size(CvrnYX,3)
    CvrnYX(:,:,iClass) = (CvrnYX(:,:,iClass)+CvrnYX(:,:,iClass)')/2;
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotGMMfit(Y,xi,XpctYX,CvrnYX,qXgivenY)

% Ns
[Nobsv,Nsamples] = size(Y);
Nclasses = size(XpctYX,2);

figure(133); clf; hold on;
switch Nobsv
    case 1
        [N,y] = hist(gather(Y),Nsamples/50);
        bar(y,N,2.2);
        pdfs = normpdf(repmat(y',[1,Nclasses]),...
            repmat(XpctYX,[length(y),1]),...
            repmat(sqrt(squeeze(CvrnYX))',[length(y),1]));
        pdfs = pdfs./sum(pdfs); % b/c (y < full support)
        pdfs = (pdfs.*xi)*sum(N);
        plot(y,pdfs);
        plot(y,sum(pdfs,2));
    case 2
        
        % colors
        if Nclasses <= 3
            colorscheme = zeros(Nsamples,3);
            colorscheme(:,1:Nclasses) = permute(qXgivenY,[2,3,1]);
        else
            colorscheme = repmat(0.5,[Nsamples,3]);
        end
        
        % plot
        scatter(gather(Y(1,:)),gather(Y(2,:)),[],colorscheme);
        scatter(gather(XpctYX(1,:)),gather(XpctYX(2,:)),'r')
        for iClass = 1:size(CvrnYX,3)
            error_ellipse(CvrnYX(:,:,iClass),XpctYX(:,iClass),...
                'conf',0.95,'style','k');
        end
    case 3
        fprintf('not yet programmed\n');
    
end
hold off;
pause(0.1);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [XNtrpYnew,CONT] = getCrossEntropy(XNtrpYold,XNtrpYnew,thr)

% cross entropy
fractionalImprovement = (XNtrpYold - XNtrpYnew)/abs(XNtrpYold);
if (fractionalImprovement < thr)
    fprintf('exiting with delta cross-entropy below tolerance\n');
    CONT = 0;
else
    CONT = 1;
end

end
%-------------------------------------------------------------------------%




















