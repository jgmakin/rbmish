%% NS219, Problem set 2
%
% This file requires the following to run:
%
% * plotKFstuff.m
%
% * plotKFstuff2.m
%
% * checkKFSteadyState.m
%
% * shadedErrorBar.m
%
% * KFsyntheticData.mat
%
% * KFdataHHS.mat
%
% Created: 05/06/13
%   by JGM



%% PROBLEM 1
% Prove that, if X and Y are independent:
%
% $$ Cov[A X + B Y] = A Cov[X] A^T + B Cov[Y] B^T,$$
%
% using the definition of covariance and of independence.



%% PROBLEM 2
% Given a dynamical system:
%
% $$ x_{t+1} = A x_t + B w_t $$
%
% $$ y_t = C x_t + r_t $$,
% 
% with known matrices $A$, $B$, and $C$; and with independent, zero-mean, 
% Gaussian white noise terms $w_t$ and $r_t$ with known covariances 
% $\Sigma_w$ and $\Sigma_r$, respectively, we are asked to solve the
% "filtering problem."
% That is, we want to know the probability (density) of the current state
% (at time $t$), given all the measurements up to and including right now:
%
% $$ p(x_t|y_0,\cdots,y_t). $$
%
% If the intial state (vector) is drawn from a Gaussian prior distribution,
% with known mean and covariance $\mu_0$ and $\Sigma_0$, respectively, then
% this posterior distribution turns out to be Gaussian, as well, so we need
% only compute its mean and covariance.  
% The celebrated computation of these cumulants, due to Kalman, takes the 
% form of an iterative update.  
% 
% In addition to the "time update" equations,
%
% $$\hat x_{t+1|t} = A \hat x_t $$
% 
% $$P_{t+1|t} = A P_{t|t} A^T + B \Sigma_w B^T, $$
% 
% the Kalman filter has a set of "measurement update" equations, which can
% be written in two different forms: the "canonical ones" that you'll see
% in most textbooks:
%
% $$ K_{t+1} = P_{t+1|t} C^T [\Sigma_r + C P_{t+1|t} C^T]^{-1} $$
%
% $$ \hat x_{t+1|t+1} = \hat x_{t+1|t} + K_{t+1}[y_{t+1} - C\hat x_{t+1|t}]
% $$
%
% $$ P_{t+1|t+1} = P_{t+1|t} - K_{t+1} C P_{t+1|t} $$
%
% and the pair we derived in class:
%
% $$ \hat x_{t+1|t+1} = P_{t+1|t+1} [C^T\Sigma_r^{-1}y_{t+1} + 
% P_{t+1|t}^{-1}\hat x_{t+1|t}] $$
%
% $$ P^{-1}_{t+1|t+1} = C^T \Sigma_r^{-1} C + P^{-1}_{t+1|t}. $$
% 
% Derive the first set of equations from the second set. 
% You will need to use the Woodbury matrix identity:
% 
% $$ (A + U C V)^{-1} = A^{-1} - A^{-1} U (C^{-1}+VA^{-1} U)^{-1}VA^{-1},$$
%
% where $A$ has size $n \times n$, $U$ is $n \times m$, $C$ is $m \times
% m$, and $V$ is $m \times n$.



%% PROBLEM 3: The Kalman Filter on synthetic data
%
% Now write a function that takes in the parameters and the emissions, and
% produces $\hat x_{t|t}$ and $P_{t|t}$ over time.  (Bear in mind that this
% stores the cumulants after the measurement update, not after the time
% update.)  Make Xhat the same size as X (i.e., Nstates x T), and make
% CvrnMat (Nstates x Nstate x T).
%
% I encoded the 

clear all;
load KFsyntheticData

%--- YOU NEED TO WRITE THIS FUNCTION ---%
% Check out the params structure: it contains all the required parameters
% for the filter---just for elegance, so we can pass fewer parameters into
% your KF function.  I also included the field 'Info0', but that's just the
% inverse of Sigma0.
KFdstrbs = KalmanFilter(params,Y);
Xhat = KFdstrbs.XHATMU;
CvrnMat = KFdstrbs.CVRNMU;
%---------------------------------------% 

% Use these to check that your filter works.

% See this m-file for detais.  Basically, the true position and covariance
% should usually fall within the posterior distribution over them.
plotKFStuff(X,Xhat,CvrnMat);

% After enough time, the filter's covariance should converge to something.
% There is a closed-form solution to this "steady-state" covariance.  The
% function below measures the difference b/n this covariance matrix and the
% one you got at time T (we assume that the simulation period is long
% enough for the filter to have reached steady state).  Thus, ressum ought
% to be very small, and if it's not, you messed up the implementation of
% the Kalman filter.
ressum = checkKFSteadyState(CvrnMat,params);




%% PROBLEM 4: The Kalman Filter on neural data (collected by Helen Shen)
%
% You've already built a Kalman Filter.  Your task is now to build a model
% of some experimental data, and then run the filter on it.  Remember, your
% model needs the following components: A, B, C, mu0, Sigma0, SigmaX, and
% SigmaY. We're going to find them one by one.
%
% The data, unfortunately, are not continuous, but broken up into trials.
% Each one consists of a center-out reach.  We're going to fit the model
% across all trials, but then test it on each trial sequentially.

clear all;
load KFdataHHS

% useful params
Nneurons = length(UnitSpikes);
Nttrials = length(S);
Ntargets = 8;
Nstates = 6;


%% The A matrix: Model the dynamics in each trial
% 
% We start by using the fact that position is just integrated velocity, and
% that velocity is just integrated position, to fill in most of the entries
% of the state-transition matrix.  That is,
%
% $$ x_{t+1} \approx x_t + \dot x_t dt $$
% 
% $$ y_{t+1} \approx y_t + \dot y_t dt, $$
%
% and so on for the velocity components.  Hence, the state transition 
% matrix is:
dt = 1/240;             % actually different for each trial, so we round
A = [1 0 dt 0 0 0;
    0 1 0 dt 0 0;
    0 0 1 0 dt 0;
    0 0 0 1 0 dt;
    0 0 0 0 0 0 ;       % These last two rows are to be determined
    0 0 0 0 0 0];       %


%% Fitting the A matrix
% To find a good A matrix, you need to use the trajectories of the arm on
% each trial.  Remember that, on average, we want 
%
% $$ x_{t+1} = A x_t. $$
%
% (We assume zero-mean noise.)
% So we collect up, for each trial, two quantities: XPast, which is all the
% dynamical information except the *last* sample of the trial; and XFuture,
% which is all the dynamical information except the *first* sample.


% init
XFuture = [];
XPast = [];
X0 = zeros(Nttrials,6);

% collect the dynamics on each trial
for iTTrial = 1:Nttrials
    thisX = [S(iTTrial).pos S(iTTrial).vel S(iTTrial).acc];
    XPast = [XPast; thisX(1:end-1,:)];
    XFuture = [XFuture; thisX(2:end,:)];
    X0(iTTrial,:) = thisX(1,:);
end



% Now find the last two rows of the A matrix---call them Aacc'

%--- YOU NEED TO WRITE THIS FUNCTION ---%
% Aacc = ...
accFuture = XFuture(:,5:6);
Aacc = (XPast'*XPast)\XPast'*accFuture;
%---------------------------------------% 

% and then fill them in:
A(5:6,:) = Aacc';
% (The transpose is a hint.)




%% Test to make sure the A matrix is any good
% This is a very weak test, but if your A matrix fails it, it is wrong.

% for a random trial....
j = ceil(158*rand);

% get the dynamics...
thisX = [S(j).pos S(j).vel S(j).acc]';

% malloc
xhat = zeros(size(thisX));

% simulate forward in time---so you're only checking your one-step
% prediction
for t = 1:(length(S(j).t)-1)
    xhat(:,t+1) = A*thisX(:,t);
end

% compare predicted and actual accelerations
figure(103); clf; hold on;
scatter(xhat(5,:),xhat(6,:))
scatter(thisX(5,:),thisX(6,:),[],'r.')
hold off;



%% p(x0): Compute the prior distribution

%--- YOU NEED TO WRITE THIS FUNCTION ---%
% mu0 = ...
% Sigma0 = ...
mu0 = mean(X0);
Sigma0 = cov(X0);
%---------------------------------------% 
% Hint: we collected above the data you'll need for this.  Remember that
% the prior is over the initial state *at the start of every trial*.


%% Bin spike counts 
% Next you're going to fit the output matrix, C.  So first you need to
% collect up spikes from the neurons into reasonable bin sizes.  Too small,
% and noise will dominate---and the assumption of Gaussian noise on the
% emissions will break.  Too large, and you'll have too few samples/trial
% to fit the C matrix.
% 
% However, using a binsize bigger than the sampling period introduces a
% complication.  If binsize = m*dt, the sampling of the dynamics will be m
% times as frequent as the sampling of the emissions (spiking rates).
% Since we used the sampling of the dynamics to fit A, we will need to
% iterate it forward m steps in time b/n emissions.  That is:
% 
% $$ x_{t+1} = A x_t,$$
% 
% so
%
% $$ x_{t+16} = A^{16} x_t.$$
% 
% So the transition matrix we use in the filter will be AA = A^(16).  We 
% shall simply average the trajectories over each of these "dark" periods.
%
% You won't need to change anything in this cell of code, but you will need
% to understand it.

% sampling rate
m = 16;                                     % 16 => 60 Hz (66.7 ms bins)
binsize = m*dt;                             %

% initialize
R = [];
Xavg = [];
XavgPast = [];                              % we'll need these later
XavgFuture = [];                            %

% loop through trials, collecting up binned spike rates
for iTTrial = 1:Nttrials
    
    % set the edges of the bins
    edges = S(iTTrial).t(1):binsize:S(iTTrial).t(end);
    Nbins = length(edges) - 1;
    
    % bin the spikes into the matrix Spks
    Spks = zeros(Nbins,Nneurons);           % malloc
    for iNeuron = 1:Nneurons
        N = histc(UnitSpikes(iNeuron).t,edges);
        Spks(:,iNeuron) = N(1:end-1);
    end
    
    % convert to spike *rate* and store
    R = [R; Spks/binsize];
    
    
    
    % To fit the C matrix, you will need both the trajectories (X) and the
    % emissions (spike rates) across each trial.  But the trajectories need
    % to be averaged across the 16-sample "dark period" when no emissions 
    % are available (a consequence of using largish bin sizes for the spike
    % counts).
    
    % get the (averaged) dynamics during this trial
    thisX = [S(iTTrial).pos S(iTTrial).vel S(iTTrial).acc]';
    thisXavg = squeeze(mean(reshape(thisX(:,1:(Nbins*m)),[Nstates,m,Nbins]),2));
    Xavg = [Xavg; thisXavg'];
    
    % You'll need these to find the transition covariance, SigmaX
    XavgPast = [XavgPast; thisXavg(:,1:end-1)'];
    XavgFuture = [XavgFuture; thisXavg(:,2:end)'];
    
    
end


% final tick
T = size(R,1);




%% The C matrix
% Now you will fit the C matrix.  Our fit will be better if we:
%
% # throw away very low firing neurons; 
%
% # take the <http://en.wikipedia.org/wiki/Variance-stabilizing_transformation square roots of the firing rates>
% (this makes them more Gaussian); and 
%
% # normalize them to zero mean and unit variance---i.e., "z-score" them.  
% Convince yourself that this doesn't lose information for the fit.
%

% massage the data
highfiringneurons = (mean(R)>1);            % greater than 5 Hz
R2 = sqrt(R(:,highfiringneurons));          % "variance stabilizing xform"
Z = zscore(R2);


% find tuned neurons
%--- YOU NEED TO WRITE THIS FUNCTION ---%
% Do something with Xavg and Z to get the C matrix
% C = ...
[beta, RsqCV, ZResCV] = linregress(Xavg,Z,'LOO');
C = beta';
figure; hist(RsqCV,20);
% Hint: when you're finished, you may need to transpose something....  
% Keep track of which dimension is where.
%---------------------------------------% 


% Find how good your fit was---I used a leave-one-out cross-validation.
% Then keep only the cells that are "tuned" (by some weak criterion:
% unfortunately, no cell will be tuned more than weakly).  My criterion
% gave me just six tuned neurons.  You can experiment with different sizes.
TunedNeurons = 1:size(Z,2); % find(RsqCV>0.15);
Y = Z(:,TunedNeurons)';
C = C(TunedNeurons,:);




%% SigmaX: Compute the transition noise
% Now we find the transition covariance matrix.  Without loss of
% generality, we can set the B matrix to the identity.  Then we need to
% find the 16-sample A matrix.

% w.l.o.g.
B = eye(Nstates);

% go m time steps before emitting
AA = A^m;


% get residuals 
%--- YOU NEED TO WRITE THIS FUNCTION ---%
% Do something with XavgPast, XavgFuture, and AA' (hint) to get the
% covariance of the transition noise.  You should also get the *mean* of
% the transition noise to make sure it's near zero for all 6 components.
% (Actually, it won't be that close for acceleration---don't worry about
% that.)
% SigmaX = ...

% get transition noise (check that the mean is ~0)
Res = XavgFuture - XavgPast*AA';
Xpinv = (XavgPast'*XavgPast)\XavgPast';
XResCV = zeros(size(Res));
for i = 1:size(XavgPast,1)
    XResCV(i,:) = Res(i,:)/(1 - XavgPast(i,:)*Xpinv(:,i));
end
% Res = XavgPast*AA' - XavgFuture;
muX = mean(Res);
SgX = cov(Res);
SigmaX = SgX;
%---------------------------------------%


%% SigmaY: Compute the emission noise
% Here you will do something very similar to the above.


%--- YOU NEED TO WRITE THIS FUNCTION ---%
% Do something with Xavg, Y', and C' (hint) to get the covariance of the
% emission noise.  You should also get the *mean* of the emission noise to
% make sure it's near zero for all components.
% SigmaY = ...
% Res = Xavg*C' - Y';
muY = mean(ZResCV);
SigmaY = cov(ZResCV);
%---------------------------------------%


%% Kalman Filter
% Now we see how well your filter does.  We really ought to test it on data
% that were held out---but then performance would be even worse, and it
% turns out that this data set doesn't give great performance.

% First we load up our parameters into the params structure, so that we
% can pass it to your Kalman Filter function

params.A = AA;
%%% params.B = B;
params.C = C;

params.mu0 = mu0';
params.Info0 = inv(Sigma0);

params.SigmaX = SigmaX;
params.SigmaY = SigmaY;

% Now we filter---separately for each trajectory, b/c otherwise the model
% will have to "jump" from target back to center, which the model can't
% possibly predict!

% init
ind = 0;
i = 1;
F = 0;
G = 0;
err = zeros(size(Y,2),Nstates);

% loop through trials
for iTTrial = 1:Nttrials
    
    % get this trial's (averaged) trajectory
    edges = S(iTTrial).t(1):binsize:S(iTTrial).t(end);
    Nbins = length(edges) - 1;
    thisX = [S(iTTrial).pos S(iTTrial).vel S(iTTrial).acc]';
    thisXavg = squeeze(mean(reshape(thisX(:,1:(Nbins*m)),[Nstates,m,Nbins]),2));
    
    % get this trial's emissions
    thisY = Y(:,(ind+1):(ind+Nbins));
    ind = ind+Nbins;
    
    % run the filter
    params.T = Nbins;   
    % (my version of the filter needs to know the final time)
    
    %--- YOU NEED TO WRITE THIS FUNCTION ---%
    % You need to insert your filter here
    KFdstrbs = KalmanFilter(params,thisY);
    Xhat = KFdstrbs.XHATMU;
    CvrnMat = KFdstrbs.CVRNMU;
    %---------------------------------------%
    
    % some checks
    [fracWithinOneStdDev,AvgVrnc] = plotKFStuff2(thisXavg,Xhat,CvrnMat); 
    i=i+1;
    F = fracWithinOneStdDev + F;        % not *quite* right, b/c the # of 
    G = AvgVrnc + G;                    % samples/trial isn't constant...
    % pause()
    
    
    % store error
    err((ind+1):(ind+Nbins),:) = (Xhat - thisXavg)';
    
end

% these should be around 0.67 (one std. dev. in each direction) 
F = F/Nttrials;

% this will ideally be "small"
G = G/Nttrials;

% find the mean square error
MSE = err'*err/(size(err,1)-1);






