% simulation of the network responses to different shifts


% init
clear; clc; close all
% load
% load results/new/wtsStdEqual
load results/numhidswts/Std050.mat
% load results/new/wtsBigSpace.mat
% load results/wts15smpllrn.mat
% load results/NN15smplLrn.mat

Nexamples = 40000;

% set gains
gain = mean([params.gmin; params.gmax]);


% test
% ErrorStats = test(wts,params,'visgain',gain(1),'propgain',gain(2));
ErrorStats = test(wts,params,'visgain',gain(1),'propgain',gain(2),...
    'propagation','Nsamples','numsamples',15);
% [ErrorStats net] = nnDecode(wts,'prop',params,'pretrained',net,...
%     'propgain',pgain,'visgain',vgain);


shiftvec = [2.5 5 7.5 10 12.5 15]; % % [5 10 15 20 25 30];
Ndims = params.Ndims;
Nmods = length(params.mods);
smin = params.smin;
smax = params.smax;
indices = reshape(1:Nmods*Ndims,Ndims,Nmods);
switch params.NS
    case 'Hand-Position'
        biasstr = 'visbias';
    case 'Joint-Angle'
        biasstr = 'propbias';
    case 'Gaze-Angle'
        biasstr = 'eyebias';
end
NSIND = find(strcmp(params.mods,params.NS));
indout = indices(:,NSIND);

% compute fisher info in whichever population is "neutral"
[D0,S0] = generateData(Nexamples,params);
SINSMerr = covInCalc(D0,S0,params);
clear D0 S0;
ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;
% *** see (1) ***

ErrorStatsArray = ErrorStats;
shiftArray = [0; 0];
nDirections = 8;
for iShift = 1:length(shiftvec)
    
    for iDirection = 0:nDirections-1;
        close all
        
        % phi = 2*pi*rand;
        phi = 2*pi*iDirection/nDirections;
        R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
        % shft = chol(ErrCovIn)'*R*[shiftvec(iShift); 0];
        shft = sqrtm(ErrCovIn)*R*[shiftvec(iShift); 0];
        
        % generate input and output data
        ErrorStats = test(wts,params,biasstr,shft,...
            'visgain',gain(1),'propgain',gain(2),...
            'propagation','Nsamples','numsamples',15);  %%% ADDED THIS LINE
%         [ErrorStats net] = nnDecode(wts,'prop',params,'pretrained',...
%             net,'propgain',pgain,'visgain',vgain,biasstr,shft);
        
        % store
        ErrorStatsArray = [ErrorStatsArray; ErrorStats];
        shiftArray = [shiftArray shft];
    end
    
end

%%%%%%%%%%
tag1 = num2str(sprintf('%02.f',gain(1)));
tag2 = num2str(sprintf('%02.f',gain(2)));
filename = ['eStats',tag1,tag2];
save(filename,'ErrorStatsArray','shiftArray','shiftvec','params');
%%%%%%%%%%


%-------------------------------------------------------------------------%
% *** (1) ***
% The input error covariance is the same for the marginal and conditional
% distributions.  Since you're using the "neutral-space population," you
% don't have to worry about transforming by Jacobians or etc.
%
% And then, you're using g = 15 to compute this cov; whereas it really
% ought to be different for the different gains (pgain = 10,20).  But then
% you couldn't compare the 30-std shift as nicely....
%-------------------------------------------------------------------------%
