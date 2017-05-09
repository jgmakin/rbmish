function KFparams = fitDynamics(S,KFparams)
% The A matrix: Model the dynamics in each trial

%-------------------------------------------------------------------------%
% Cribbed: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%

% params
Nstates = KFparams.Nstates;
Ndims = KFparams.Ndims;


% initialize the A matrix
dt = KFparams.dt;
A = diag(ones(Nstates,1)) + diag(dt*ones(Nstates-Ndims,1),Ndims);
A(end-Ndims+1:end,:) = 0;   % last Ndims rows are t.b.d


% initialize
Nttrials = length(S);
endinds = zeros(Nttrials,1);
ind = 0;
for iTTrial = 1:Nttrials
    Nsamples = length(S(iTTrial).t);
    ind = ind+Nsamples;
    endinds(iTTrial+1) = ind;
end
X = zeros(ind,Nstates);



% collect the dynamics on each trial
for iTTrial = 1:Nttrials
    switch Nstates/Ndims %%% not very elegant
        case 1
            thisX = [S(iTTrial).pos];
        case 2
            thisX = [S(iTTrial).pos S(iTTrial).vel];
        case 3
            thisX = [S(iTTrial).pos S(iTTrial).vel S(iTTrial).acc];
    end
    X((endinds(iTTrial)+1):(endinds(iTTrial+1)),:) = thisX;
end
XFuture = X(~ismember(1:size(X,1),endinds+1),:);
XPast = X(~ismember(1:size(X,1),endinds),:);
X0 = X(endinds(1:end-1)+1,:);


% Now find the last Ndims rows of the A matrix---call them Aacc'
accFuture = XFuture(:,end-Ndims+1:end);
[Aacc, RsqCV] = linregress(XPast,accFuture,'LOO');
fprintf('row %i fit with R^2 = %0.3f\n\n',[1+Nstates-(Ndims:-1:1);RsqCV]);
A(end-Ndims+1:end,:) = Aacc';

% get transition noise (check that the mean is ~0)
Res = XFuture - XPast*A';
Xpinv = (XPast'*XPast)\XPast';
ResCV = zeros(size(Res));
for i = 1:size(XPast,1)
    ResCV(i,:) = Res(i,:)/(1 - XPast(i,:)*Xpinv(:,i));
end


% noise terms
muX = mean(ResCV);
SigmaX = cov(ResCV);

% p(x0): Compute the prior distribution
mu0 = mean(X0);
Sigma0 = cov(X0);



% load into params
KFparams.A = A;
KFparams.SigmaX = SigmaX;
KFparams.Sigma0 = Sigma0;
KFparams.Info0 = inv(Sigma0);
KFparams.mu0 = mu0';
KFparams.muX = muX';


end





