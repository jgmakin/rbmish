function [statsL,statsN] = condInfCondErrs(wts,params)
% The results of "conditional inference"---inferring VIS from PROP and vice
% versa---can't be averaged over the whole space, b/c the different priors
% in the different spaces bias the answers in different directions.  That's
% fine: this function just shows conditional errors for a few locations
% (see the data structure proplocs).
%
% For simplicity you also keep the gain fixed at each location
% (params.swing = 0).
%
% Empirical investigation reveals that a single updown pass suffices to
% generate a good estimate of the hand's visual coordinates from its
% proprioceptive ones (and vice versa), so this function doesn't bother to
% clamp and iterate; for that, see condInfPP.m
%
% USAGE:
%   load results/numhidswts/Std050
%   stats = condInfOnePass(wts,params)
%
% NB YOU SHOULD PROBABLY FIRST COMMENT OUT THE FOLLOWING:
%     "center of mass detects a zero" in CoM.m


%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S etc.
% Revised: 12/10/13
%   -changed indexing of stats, based on new output format from
%   estStatsCorePP.m
% Created: 10/19/12
%   by JGM
%-------------------------------------------------------------------------%

%%%%%% TO DO %%%%%%
% (1) fix (esp. plotting?) to work with length(params.mods) = 3, so that you can check
%   coordinate transformation
%%%%%%%%%%%%%%%%%%%

% init
gains = mean([params.gmin; params.gmax]);
params.gmin = gains;
params.gmax = gains;
Nexamples = 2000;

% set the locations to be investigated (as a fraction of prop space)
% proplocs = [0.5 0.5; 0.25 0.25; 0.25 0.75; 0.75 0.25; 0.75 0.75]';
proplocs = [0.1 0.1; 0.9 0.9];

% get the conditional errors for conditional inference
[statsL,statsN] = getCondInfCondErrStats(proplocs,Nexamples,wts,params);

% plot
plotCondInfStats(statsN,proplocs,Nexamples,params);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [statsI,statsF] =... % [statsL statsN] =...
    getCondInfCondErrStats(proplocs,Nexamples,wts,params)

% init
HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED, pool = parpool(4); end        % open pool once for all
numlocs = size(proplocs,1);                     % how many locs do we test?
Nmods = length(params.mods);

% loop through the locations!
for loc = 1:numlocs
    
    % transform this fraction of prop space into a set of joint angles
    p0.mu = scalefxn(proplocs(loc,:),[0;0],[1;1],...
        params.roboparams.thmin,params.roboparams.thmax);
    p0.cov = 0;
    
    % generate data at one point
    [D0,S0] = generateData(Nexamples,params,'stimulusprior',p0);
    
    % sequentially turn off modalities (cond. inf. in both directions)
    for modeOFF = 1:Nmods
        
        [statsI{loc}{modeOFF},statsF{loc}{modeOFF}] =...
            condInfPP(D0,S0,wts,params.mods{modeOFF},params);
        
%         % zero out one input modality, for all examples
%         [Di S.zvec] = killoneinput(D0,params.mods{modeOFF},length(params.mods));
%         
%         % stats for a single up-down pass (use neutral-space stats only)
%         [~, Do] = updown(Di,wts,params,'means','quiet');
%         [statsL{loc}{modeOFF}, statsN{loc}{modeOFF}] =...
%             estStatsCorePP(x0,params,'CoM',Di,Do);

    end
    
end
if HADBEENCLOSED, delete(pool); end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotCondInfStats(stats,proplocs,M,params)

% init
numlocs = size(proplocs,1);
Nmods = length(params.mods);
thmin = params.roboparams.thmin;
thmax = params.roboparams.thmax;
setColors;

% dimensionality of the stimulus?
switch params.Ndims
    case 1  % 1D vars
        %%% this case may still be broken and want fixing
        
        % inputvar = statsOnePass{1}.cov + statsOnePass{3}.cov;
        % outputvar = statsFinal{length(params.mods) + 2}.cov;
        
        for loc = 1:numlocs
            for modeOFF = 1:Nmods
                othermods = setdiff(1:Nmods,modeOFF);
                
                cov = 0;
                for j = 1:length(othermods)
                    cov = cov + stats{loc}{modeOFF}{othermods(j),1}.cov;
                end
                cov
                stats{modeOFF}{modeOFF,2}.cov
                
                pause;
                fprintf('\n');
            end
        end

    case 2 % 2D vars
        
        % place the figure
        close all;
        h = figure;
        set(h,'position',[213 261 1028 450]);
        titlestr = {'forward kinematics','inverse kinematics'};
        
        for loc = 1:numlocs
            
            th = scalefxn(proplocs(loc,:),[0;0],[1;1],thmin,thmax);
            for modeOFF = 1:Nmods
                modeON = mod(modeOFF,2)+1;
                subplot(1,2,modeOFF);
                hold on
                plotErrStats(stats{loc}{modeOFF}{modeON,1},M,OPTcolor,th);
                plotErrStats(stats{loc}{modeOFF}{modeOFF,2},M,EFHcolor,th);
                title(titlestr{modeOFF});
                xlabel('$\theta_1$','Interpreter','latex');
                ylabel('$\theta_2$','Interpreter','latex');
                legend('input','inferred');
                axis equal;
                hold off
            end
            
        end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotErrStats(stats,M,clr,th)
%%% change color selection to use RBMcolor, etc. etc.

cntr = stats.mu + th;
cvr = stats.cov;
h = error_ellipse(cvr,cntr,'conf',.95); 
set(h,'LineWidth',1,'Color',clr);
h = error_ellipse(cvr/M,cntr,'conf',.95);   
set(h,'LineWidth',1,'Color',clr);

end
%-------------------------------------------------------------------------%