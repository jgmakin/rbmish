% KF4HHS
%   A Kalman filter for HHS's data.  You need to distribute with this the
%   following files:
%       plotKFstuff.m
%       plotKFstuff2.m
%       checkKFSteadyState.m
%       shadedErrorBar

%-------------------------------------------------------------------------%
% Revised: 11/25/13
%   -fixed yet another bug: you were fitting by minimizing MSE (that's what
%   linrgnLOOCV.m does), but not using a column of ones!!  Now you're
%   demeaning both X and Y before fitting C (and the Xval MSE is lower).
% Revised: 05/20/13
%   -cleaned up; eliminated switch statement and multiple versions of files
%   (binSpikeCounts1, 2, etc.) in favor of switch statements *within* those
%   files.
%   -corrected a few bugs in the m-files called by this fxn
% Revised: 05/17/13
%   -changed to use one of three different methods: nonoverlapping windows,
%   sliding windows, or downsampled data.
% Revised: 05/13/13
%   -improved, functionized
% Revised: 05/06/13
% Revised: 05/03/13
% Revised: 05/02/13
% Revised: 05/01/13
% Created: 04/30/13
%   by JGM
%-------------------------------------------------------------------------%



% load and init
clear all;
datadir = 'C:\#code\HHS\extracteddata\';
% tag = 'D080616';
tag = 'D080620';
load([datadir,'KFtuningdataHHS',tag]); % "the tuning series"

% useful params
KFparams.Nstates = 6;
KFparams.Ndims = 2;
KFparams.m = 16;                              % 16 => 60 Hz (66.7 ms bins)
KFparams.dt = 1/240;                          %

% get all the parameters
KFparams = fitDynamics(St,KFparams);
testAmatrix(KFparams.A,St);

KFparams.BINMETHOD = 'nonoverlappingwindow';
% 5.1742e+04, 5.8648e+04                'slidingwindow'
% 2.4413e+04, 3.7337e+04                'nonoverlappingwindow'
% 2.5023e+04, 2.6811e+04                ", corrected (<emis. noise>=0)
% 2.0733e+04, 2.5958e+04                ", corrected (used col of ones)
% 2.0963e+04, 2.5243e+04                ", fit whole A to downsampled data
% 3.1441e+04, 4.7949e+04                'downsampled'


% fit a fully observed LDS
[R,Xout,endinds] = binSpikeCounts(St,UnitSpikesT,KFparams);
KFparams = refitTransitionNoise(Xout,endinds,KFparams); 
vec = mean(R)>0;
KFparams = fitEmissions(Xout,R(:,vec),KFparams);
Y = R(:,vec);
        
% run the KF
[MSE, Rsq] = filterAllTrials(Xout,Y,endinds,KFparams);
det(MSE(1:2,1:2))



if 1
    % run this Kalman filter on another data set ("the target series")
    load([datadir,'KFreachingdataHHS',tag]); 
    [Rr,Xoutr,endindsr] = binSpikeCounts(Sr,UnitSpikesR,KFparams);
    %%% Yr = Rr(:,vec) - repmat(params.muYX,size(Rr,1),1);
    Yr = Rr(:,vec);
    clear Rr
    [MSEr, Rsqr] = filterAllTrials(Xoutr,Yr,endindsr,KFparams);
    
    det(MSEr(1:2,1:2))
end






%%%%%%%
load(['KFdataJEO','I121204']);
clear paramsJEO;
paramsJEO.Nstates = 6;
paramsJEO.m = 16;                              % 16 => 60 Hz (66.7 ms bins)
paramsJEO.dt = 1/240;  
paramsJEO.BINMETHOD = 'nonoverlappingwindow';
paramsJEO = fitDynamics(S,paramsJEO);
[Rjeo,Xoutjeo,endindsjeo] = binSpikeCounts(S,UnitSpikes,paramsJEO);
paramsJEO = refitTransitionNoise(Xoutjeo,endindsjeo,paramsJEO);
vecJEO = mean(Rjeo)>0;
paramsJEO = fitEmissions(Xoutjeo,Rjeo(:,vecJEO),paramsJEO);
Yjeo = (Rjeo(:,vecJEO) - repmat(paramsJEO.muYX,size(Rjeo,1),1))';
[MSEjeo, RsqJEO] = filterAllTrials(Xoutjeo,Yjeo,endindsjeo,paramsJEO);
RsqJEO
%%%%%%%%


% shuffle the indices of R
[foo,a] = sort(rand(1,size(Rjeo,1)));
[foo,b] = sort(rand(1,size(Rjeo,2)));
RR = Rjeo(a,b);
Rjeo = RR;
%%%%%%%%




%%%
% (1) try building C matrix on averaged data, but testing on windowed data??

%%%
%%% downsampled (leave out last) 2.4875e+04
%%% downsampled (leave out first)  2.0848e+04




% best: 2.4278e+04 (m=16, no sqrt transform, mean(R)>0)
% Rsq' = [0.8914    0.9106    0.6173    0.7177    0.2090    0.3410]
    
    
%%% det:  (m=16)
%%% det: (m=24)
%%% det: (m=30)
%%% det: (m=36)

%%% --(1) center X at zero??
%%% --(2) include the extra ones in X??
%%% --(3) PCA the firing rates??
%%% --(4) non-zero-mean noise??
%%% --(5) only fit 2x2 A matrix (acc -> acc)??










