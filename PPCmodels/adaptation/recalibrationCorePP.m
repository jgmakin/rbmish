function [IntegL0,DoBF,stVWDR,stWDR,stEMP] =...
    recalibrationCorePP(Nsubjects,Nbatches,shft,wts,params)
% the inner recalibration loop (over subjects); gets called by
% recalibration3 and recalibration4.

%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed the indexing of IntegL0 based on new indexing of shatL in
%   estStatsCorePP.m
% Cribbed: 10/02/12
%   -from recalibration3.m
%   by JGM
%-------------------------------------------------------------------------%


% keep things clean
close all; clc;

% generate biased data
[Di,Si] = generatebiaseddata(Nbatches,shft,params);

% do the first pass to find where the original integrated estimates are
IntegL0 = getInitialIntegs(Di,Si,wts,params);

% recalibrate
fprintf('\n\nAveraging across %i subjects...\n\n',Nsubjects);
DoBF = 0;
stEMP = 0;

% ready the parallel processors
[pool,HADBEENCLOSED] = parallelInit;
spmd
    RandStream.setGlobalStream ...
        (RandStream('mt19937ar','seed',sum(100*(1 + labindex/numlabs)*clock)));
end

% loop through all subjects
parfor i = 1:Nsubjects
    [thisDoBF,stepsEMPAvg] = recalibrate(Di,IntegL0,wts,params);
    DoBF = thisDoBF + DoBF;
    stEMP = stepsEMPAvg + stEMP;
end
if HADBEENCLOSED, delete(pool); end
DoBF = DoBF/Nsubjects;
stEMP = stEMP/Nsubjects;


% plot the adaptations
% [stVWDR,stWDR] = adaptByRule(Di,reshape(xi(1,:,1),2,2),50,0.01,1000,params);
[stVWDR,stWDR] = adaptByRule(Di,squeeze(Si(1,:,:)),50,0.01,1000,params);
%%% or 100, 1, 1000.  Might want to make these settable by a varargin.
plotAdaptations(Si,IntegL0(:,:,:,end),DoBF,stVWDR,stWDR,stEMP,params);



end