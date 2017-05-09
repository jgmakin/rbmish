function params = fitEmissions(X,Y,params)
% Fit (leave-one-out cross-validated linear regression) the emission matrix
% (C), the emission noise (SigmaY), and the emission-noise mean (i.e. the
% residual mean, muY).
%
% NB that this assumes the input argument X has not been augmented with a 
% column of ones, and adds it internally.

%-------------------------------------------------------------------------%
% Revised: 11/25/13
%   -added the column of ones to X (don't worry! this just pulls out the
%   mean of the emissions, which you subtract off in KF4HHS.m (b/c your KF 
%   assumes zero-mean emissions).
% Cribbed: 05/16/13
%   -from KF4HHS.m
%   by JGM
%-------------------------------------------------------------------------%


% fit
[beta, RsqCV, ResCV] = linregress([X ones(size(X,1),1)],Y,'LOO');

% plot
figure(472); clf; hist(RsqCV,20);


% you may want to pare down Y...
C = beta(1:end-1,:)';

% SigmaY: Compute the emission noise
switch params.BINMETHOD
    case 'slidingwindow', SigmaYX = cov(ResCV)*params.m; %%% a hack, de obvio
    otherwise, SigmaYX = cov(ResCV);
end

% store
params.C = C;
params.muYX = beta(end,:)';
params.SigmaYX = SigmaYX;

end