function [MSE, Rsq] = filterAllTrials(X,Y,endinds,params)
% run the filter on each trial separately

%-------------------------------------------------------------------------%
% Revised: 11/25/13
%   -changed Y to Y' (now matches X)
% Cribbed: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%



% init
Nttrials = length(endinds)-1;
Nstates = size(params.A,1);
Nsamples = endinds(end);
i = 1; F = 0; G = 0;

% malloc
err = zeros(Nsamples,Nstates);


% loop through trials
for iTTrial = 1:Nttrials
      
    % get this trial's data
    thisY = Y((endinds(iTTrial)+1):(endinds(iTTrial+1)),:)';
    thisX = X((endinds(iTTrial)+1):(endinds(iTTrial+1)),:)';
   
    % run the filter
    params.T = endinds(iTTrial+1)-endinds(iTTrial);
    KFdstrbs = KalmanFilter(params,thisY);
    Xhat = KFdstrbs.XHATMU;
    CvrnMat = KFdstrbs.CVRNMU;
    
    % some checks
    if 0
        [fracWithinOneStdDev,AvgVrnc] = plotKFStuff2(thisX,Xhat,CvrnMat);
        i=i+1;
        F = fracWithinOneStdDev + F;        % not *quite* right, b/c the # of
        G = AvgVrnc + G;                    % samples/trial isn't constant...
        pause();
    end
    
    % store error
    err((endinds(iTTrial)+1):(endinds(iTTrial+1)),:) = (Xhat - thisX)';
    
end

% these should be around 0.67 (one std. dev. in each direction) 
F = F/Nttrials;

% this will ideally be "small"
G = G/Nttrials;

% find the mean square error
MSE = err'*err/(Nsamples-1);
MST = X'*X/(Nsamples-1);
Rsq = 1 - diag(MSE)./diag(MST);


end