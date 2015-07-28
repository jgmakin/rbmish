% LDStest   Test filters for linear dynamical systems
%   A script

%-------------------------------------------------------------------------%
% Revised: 01/29/15
%   -pushed looping over "experiments" and model types ("filters") entirely
%   into multiXprmtStats.m.
% Revised: 01/28/15
%   -completely revised, cleaned up, to run with many experiments and error
%   bars.
% Created: 09/08/14
%   by JGM
%-------------------------------------------------------------------------%


clear; clc;

tic; 

% init
params = setParams;

% different results for different models
switch params.MODEL
    case {'1DrEFH'} % '2DrEFH','1DtRBM'}
        filterNames = {'sensory','EM1stOrd','rEFH','EM2ndOrd','KFtrue'};
    case '1DrEFHwithEC'
        filterNames = {'sensory','KFobsNoCtrl','EM2ndOrd','rEFH','EM3rdOrd','KFtrue'};
        % filterNames = {'EM2ndOrd','yrtest'};
end
Nfilters = length(filterNames);
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;
%%%params.Ncases = 320;

[MedianMSEs,MSEerrorBars,xaxislabels,clrNames] =...
    multiXprmtStats(filterNames,params);

% plot
for iMod = 1:params.Nmods
    nonNanInds = ~isnan(MedianMSEs(:,iMod));
    tikzBarGraph(MedianMSEs(nonNanInds,iMod),...
        MSEerrorBars(nonNanInds,:,iMod),...
        2,0.32,xaxislabels(nonNanInds,iMod),'estimator',...
        'mean square error','',clrNames(nonNanInds,iMod),...
        ['1DerrorStats',params.mods{iMod},date]);
    
end
        

toc;
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
if 0
    
t0 = 1;

% errors and estimates vs. time; th1 vs. th2
plotFilterErr(t0,LDSdataTest,params,...
    pSENSORY,pKFEM,pKFobsNoU,pKFobs,pKFtrue,pEFH);
% plotFilterErr(t0,1:size(S0,1),S0,params,pSENSORY,pKF);

% plot the arm being tracked
j = ceil(rand*params.Ncases);
plotArmTracking(j,S0,params,pSENSORY,pRBM,pKF);

end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
if 0
% plot center of mass of putative messed up ones (due to wrapping)
Ndims = params.Ndims;
N = params.N;

for i=1:1000
    foo0 = reshape(R0(j,:,i)',[N*ones(1,Ndims),1]);
    foo1 = reshape(R1(j,:,i)',[N*ones(1,Ndims),1]);
    figure(1); imagesc(foo0); title(num2str(i));
    figure(2); imagesc(foo1); title(num2str(i));
    pause()
end

end
%-------------------------------------------------------------------------%




%

% TO DO

% -(1) make the trajectories bounce?
% -(1.5) train a network on that...
% -(2) make a "movie" of the arm being tracked
% -(3) "fix" the decoding of the output (center of mass on a torus)
% -(4) write the optimal tracking in the case of bouncing...
% -(5) Kill the "new" recurrency, just use the downpass....
%   [doesn't work so hot]
% -(6) think about dialing down the "final momentum," which induces a huge
% error every time it kicks in.
%   [You've reworked the way momentum is calculated, essentially making the
%   gradient pass through a low-pass filter!]
% -(7) Neither bouncing nor wrapping works that well: both introduce errors
% into the RBM.  So you've written dynamics in which the walls repel.  
%   [But this, altho' better on avg than the input, is also suboptimal.]
%   [The badness appears to accrue on trajectories (cases) that near the
%   edges of the workspace.  How now?]
% (8) try pKF2 ????
% -(9) Tried recursing with V1 rather than V0---or really, V0bar.  And lo,
% it learns just as well----except that in *testing*, one must recurse with
% V0bar to get the best results (even tho' it was trained on V1!!!).
% -(10) Try a differently nonlinearity: a quadratic bowl (centered at the
% center of the workspace)
%   [Works well, though not obviously more interesting than the repelling
%   walls.  In particular, both the KF and the dRBM do slightly better.  On
%   the other hand, the one big advantage is that it doesn't ever leave the
%   feasible area and have to get reset.  This should really speed up
%   training, and it makes the distributions more sensible (the resets are
%   essentially performing some kind of rejection sampling...).
% -(11) Compared with the KF learned via EM rather than guessed.  This KF
% is much better :-(
% -(12) Increase number of units??
%   [Went up to [1575 x 1350], didn't really help.  Will try even more.]
% -(13) Tried using V0 rather than V0bar in the recusion for training.
%   This seems to work just as well.
% -(14) I think you want to initialize the dRBM differently during testing.
%   Specficially, don't use the zeros vector, which is fucked up; do
%   conditional inference or something first.... [see (20)]
% (15) Use actual reach trajectories, say from JEO's data!!
% (16) Back to applying the dRBM to real data.  I think you should really
%   try using regression to decode the hidden units, rather than the 
%   RBM-type method.  At least at first....
% -(17) PNS's idea: use a control to steer the system away from the
%   boundaries, eliminating boundary issues.  Then the dRBM has to learn to
%   use the control, as well....
% -(18) Related: update yr getLDSparams and etc. to account for KF with
%   input!!
% ?(19) Related: try putting (17) into 1D, b/c it takes too long to train
%   this network; or more precisely: to determine if extra units are 
%   helping (b/c with 1e7 weights, the network probably requires more 
%   training....)
% -(20) Eliminate training of the V(t-1) to V(t) weights on the first step
%   of every trajectory (see your proof).
%  (21) Write a little code to *premake the training data*, and then just
%   load it in RBM rather than recreate it!
% ?(22) issue: for the EM-trained LDS, a version with an input does no 
%   better than one with no input!!  And worse than one with an input 
%   learned with fully observed data.
%   [see working notes (216)]
% (23) can you compute (fake) cross-entropies for the EFH?  You would need
%   to infer the cumulants of the marginal distribution of---the center of
%   mass.



% (26) re-run EM-trained model with data with higher gain (yielding
% reliabilities like those for the rEFH)
% -(27) consider replacing all Shat with Xpct (and etc.)
% (28) Consider replacing 40 trajectories of 1000 time steps apiece with 1
% trajectory of 40,000 time steps.  This could eliminate the weirdness of
% starting with a vector of all zeros/blank input for the first input.  B/c
% the minibatches are trained sequentially, you could just seed each one
% with the last one's ending vector!  ...And then you could KF on the whole
% trajectory of 40,000 time steps.
%   If you're right that the model is trying to remember all of its inputs
%   for all time, this would test that--since it presumably can't do that
%   for 40,000 time steps w/only 225 feedback units.  Of course, it can't
%   do that for 1000 time steps either, but then it isn't perfect....

% -(29) Massive change: re-implement "noisily observed control" as another
% state, thus making the change not to the filter and smoother, but to the
% updater, e.g. right after the E step.
% -(30) See if you can't find a way around requiring params.swing = 0 for
% wrapped data
% [done]
% -(31) EM4LDS.m probably shouldn't quit on the first epoch for a given data
% set.
% (32) see if you can explicitly store the data (Poisson, Bern) as ints??
% (33) funtionize, a la mastertest.m?







 % params.dynamics.meta = 'RandInitWithEC';
 % LDSparamsEMNoCtrlDynamics = EM4LDS(size(params.dynamics.G,1),params);
 % load('dynamical\LDSparamsEM1DrEFHwithECRandInitWithEC.mat') %140904
 % load('dynamical\141003\LDSparamsEM1DrEFHwithECRandInit2ndOrder.mat') %141003
 %%% load('dynamical\141003\LDSparamsEM1DrEFHwithECRandInit2ndOrder141007.mat')
 %%% LDSparamsEMNoCtrlDynamics = LDSparamsEM;
 % params.dynamics.meta = 'RandInit';
 % LDSparamsEM = EM4LDS(size(LDSdata.Z,2),params);
 % load('dynamical\LDSparamsEM1DrEFHwithECRandInit.mat')       %140904
 % load('dynamical\141003\LDSparamsEM1DrEFHwithECRandInit.mat')       %141003
 %%% load('dynamical\141003\LDSparamsEM1DrEFHwithECRandInit141007.mat')
 
 
 
 
 
 
 
 
 
 % load a trained RBM?
% load('results/finalwts/wts1Dcontrolled140428','wts','params');
% load('results/finalwts/wts1Dwrapping140519','wts','params');
% load('results/finalwts/wts1dControlled140428.mat','wts','params');
% load('results/finalwts/1DwrappingHVNdamped.mat')
% load('results/nonfinalwts/wts1Dcontrolled140904.mat')
% load('results/nonfinalwts/wts1Dcontrolled140918.mat')
% load('results/nonfinalwts/wts1Dcontrolled141006.mat')


% Final results (10/14, 1/15)
% load('dynamical/finalwts/wts1DrEFH141029.mat')
% load('dynamical/finalwts/wts1DrEFHwithEC141027.mat');
% load('dynamical/finalwts/wts1DrEFHManyRuns150120.mat')
% load('dynamical/finalwts/wts1DtRBM150116.mat')
%%% load('results\nonfinalwts\LDSparamsEM1DrEFHwithECRandinit3rdOrd141028.mat')
% params.Ncases = 100;

% load('C:\Users\makin\Desktop\caca\wts1DrEFHoverdamped141030.mat')
% load('C:\Users\makin\Desktop\caca\LDSparamsEM1DrEFHRandInit1stOrdoverdamped.mat')
% LDSparamsEM1stOrd = LDSparamsEM;
% load('C:\Users\makin\Desktop\caca\LDSparamsEM1DrEFHRandInit2ndOrdoverdamped.mat')
% LDSparamsEM2ndOrd = LDSparamsEM;


%%% load('results/nonfinalwts/wts1D3rdOrd7.mat');
% load('results/nonfinalwts/wts1D3rdOrd7XXX.mat');
% load('C:\Users\makin\Desktop\caca\wts1D3rdOrd5XXX.mat');
% load('results\nonfinalwts\wts1D3rdOrd141014.mat');
% load('C:\Users\makin\Desktop\caca\wts1D3rdOrd141014b.mat');
% load('C:\Users\makin\Desktop\caca\wts1D3rdOrd141015b.mat');

% load('C:\Users\makin\Desktop\caca\wts1D3rdOrd141017.mat');