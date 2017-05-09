% STANDARDCOND2
%   STANDARDCOND2 produces (from scratch, or from a file, depending on the
%   RESTART param,) the conditional statistics in the standard model, at
%   stimuli that uniformly (grid-wise) cover the space of angles.  (These
%   results essentially subsume those in STANDARDCOND.m.)  Then it plots
%   them in some nice way.



clear; clc; close all
GAUSSIANPRIOR = 0;
N = 15;

if ~GAUSSIANPRIOR
    % wtsfile = 'wtsStdEqual'; % 'results/new/wtsBigSpace.mat'; % 'results/wts15smpllrn.mat';
    % wtsfile = 'results/numhidswts/Std050.mat';
    wtsfile = 'results/091012/wtsStd.mat';
    % wtsfile = 'results/new/wtsAdd90200';
    condfile = 'results/new/condStats.mat';
    plotvec = [4*N+(N+1)/2 N*(N-1)/2+2 N*(N-1)/2+(N-1) 10*N+(N+1)/2];
else
    wtsfile = 'results/new/wtsPrior.mat';
    condfile = 'results/new/condStatsPrior.mat';
    plotvec = [6*15+8 7*15+7 7*15+9 8*15+8];
end
load(wtsfile);
RESTART = 1;
NSIND = 2;

%%
if RESTART
    k = 1;
    ErrorStatsArray = [];
%     for i = linspace(params.eyemin,params.eyemax,N)
    for i = linspace(params.roboparams.thmin(1),params.roboparams.thmax(1),N)
        for j = linspace(params.roboparams.thmin(2),params.roboparams.thmax(2),N)
            % for j = linspace(params.thmin(2)+0.07,params.thmax(2)-0.07,N)
            
            p0.mu = [i; j]; % i
            p0.cov = zeros(params.Ndims);
            
            ErrorStats = test(wts,params,'stimulusprior',p0,...
                'propagation','Nsamples','numsamples',15);
            % ErrorStats = test(wts,params,'prior',p0);
            close all;
            ErrorStatsArray = [ErrorStatsArray; ErrorStats];
            
            % caca(:,k) = p0.mu;
            k=k+1;
        end
    end
    filename = ['conds',date];
    save(filename,'ErrorStatsArray','params');
else
    
    load(condfile);
end

%%
% close all;
figure(1); hold on;
k = 0;
setColors
decodeColor = EFHcolor;
clr = [decodeColor; OPTcolor];      
Nexamples = 40000;

for i = linspace(params.roboparams.thmin(1),params.roboparams.thmax(1),N)
    for j = linspace(params.roboparams.thmin(2),params.roboparams.thmax(2),N)
        % for j = linspace(params.thmin(2)+0.07,params.thmax(2)-0.07,N)
        k=k+1;
        locP(:,k) = [i; j];
        locV(:,k) = FK2link([i; j],params.roboparams,1);
        MargErrors = ErrorStatsArray(k,[1,3:4]);
        for q = 3:-1:2 % 4:-1:3
            e = MargErrors{q};
            
            covP = e{NSIND}.cov;
            biasP(:,k) = e{NSIND}.mu;
            % covV = e{1}.cov;
            biasV(:,k) = e{1}.mu;
            
            figure(1);
            h = error_ellipse(covP,biasP(:,k)+[i;j],'style','b','conf',0.95);
            set(h,'LineWidth',1,'Color',clr(q-1,:)); %,'LineStyle',estilo);
            
        end

        if sum(k == plotvec)
            figure(1);
            h = error_ellipse(covP,biasP(:,k)+[i;j],'style','r','conf',0.95);
            dispErrCovs(MargErrors,Nexamples,params);
            %%% 1000 is hard-coded to match the usual numbers in test....
        end
    end
end
hold off;
figure
quiver(locP(1,:),locP(2,:),biasP(1,:),biasP(2,:),1)
figure;
quiver(locV(1,:),locV(2,:),biasV(1,:),biasV(2,:),1)


% 







