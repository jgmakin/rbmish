function [UKFparams, XhatUKF, RsqUKF] = UKFDecoder(Rtrain,Xtrain,T,...
    Rtest,Xtest,SStot,varargin)
% function [UKFparams, RsqUKF] = UKFDecoder(Rtrain,Xtrain,T,Rtest,Xtest,...
%   SStot,Ntapspast,Ntapspast);
%
% UKFdecoder wrapper file
%
% Rtrain - training spike rates  - Nsamples x Nneurons
% Xtrain - traning kinematics    - Nsamples x Nstates
% T      - num training samples  - Nsamples x 1
% Rtest  - testing spike rates   - Nsamples x Nneurons
% Xtest  - testing kinematics    - Nsamples x Nstates
% SStot  - sumsquares of testing - 1        x Nstates
%          SStot = sum((Xtest - mean(Xtest)).^2);
% tvec   - unused

%-------------------------------------------------------------------------%
% Revised: 04/18/17
%   -added arguments, code to control the tap order
% Created: 04/13/17
%   by JEO
%-------------------------------------------------------------------------%

% taps, etc.
fparams = [];  % extra params for f, not needed
ftaps = defaulter('ftaps',2,varargin{:});       % future taps
ptaps = defaulter('ptaps',3,varargin{:});       % past taps
mtaps = defaulter('mtaps',1,varargin{:});       % backbone; <=ftaps+ptas
htaps = defaulter('htaps',1,varargin{:});       % spiking history taps
alphamax = 16; % maximum regularization (ridge regression) parameter

% f: backbone (px py vx vy ax ay) -> variables encoded linearly by spikes
f = @(a,b)([...
    a;...
    sqrt(a(1,:).^2+a(2,:).^2);...
    sqrt(a(3,:).^2+a(4,:).^2);...
    sqrt(a(5,:).^2+a(6,:).^2);...
    a(1,:).*a(3,:);...
    a(2,:).*a(4,:);...
    ]);

% fit
params = fit_ar_ukf_hist_rrcv(gather(Xtrain),gather(Rtrain),f,fparams,...
    ftaps,ptaps,mtaps,htaps,alphamax);

% re-fit initial states from *all* training data
xinit = gather(mean(Xtrain));
vinit = gather(cov(Xtrain));

% test
[XhatUKF, ~] = ar_ukf_hist(gather(Rtest),params,xinit,vinit,1);

% compute R^2 using predictions for *current* kinematic vars
Nstates = size(Xtrain,2);
indsNow = (1:Nstates) + (ptaps-1)*Nstates;
RsqUKF = 1 - sum((XhatUKF(:,indsNow)- Xtest).^2)./SStot;

% todo fill these
UKFparams.A = params.F;
UKFparams.SigmaX = [];
UKFparams.nuX = [];
UKFparams.C = [];
UKFparams.SigmaYX = [];
UKFparams.muYX = [];
UKFparams.mu0 = [];
UKFparams.Cvrn0 = [];
UKFparams.T = [];
UKFparams.params = params;

end
