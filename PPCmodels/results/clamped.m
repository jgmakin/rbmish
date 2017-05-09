%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed indexing of statsOne and statsInf based on new output format
%   from estStatsCorePP.m.
% Revised: 08/22/12
%   -fixed results display for 1D inputs
% Revised: 08/16/12
%   -
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

clear; clc; % close all

% init
RESTART = 1;
filename = 'clampedstats';

% filename = 'clampedstatsStd';
% load results/numhidswts/Std050.mat;
% modstrs = {'vis','prop'};

% filename = 'clampedstatsAdd';
% load results/addwts/add180160.mat
% modstrs = {'vis','prop','eye'};

Nmods = length(params.mods);

% generate data
%%%%%
% p0.mu = (params.thmax + params.thmin)/2;
% p0.mu = [0;1];
% p0.cov = zeros(params.Ndims);
% [D0 x0] = generateData(10,params,'visgain',gain(1),'propgain',gain(2),...
%     'stimulusprior',p0);
% %%%%%
% [D0 x0] = generateData(10,params,'visgain',gain(1),'propgain',gain(2));
Nexamples = 20000;
[D0,S0] = generateData(Nexamples,params);

% get stats
if RESTART
    
    % conditionally infer, first in one direction then the other
    for modeOFF = 1:Nmods
        tic
        [statsOne{modeOFF},statsInf{modeOFF}] =...
            condInfPP(D0,S0,wts,modstrs{modeOFF},params);
        toc
    end
    save(filename,'statsOne','statsInf')
    
else
    load(['results/new/',filename])
end




%%

switch params.Ndims
    case 2
        % place the figure
        close all;
        h = figure;
        set(h,'position',[213 261 1028 450]);
        titlestr = {'forward kinematics','inverse kinematics'};
        Nexamples = size(D0,1)*size(D0,3);
        
        % plot
        for modeOFF = 1:Nmods
            
            modeON = mod(modeOFF,2)+1;
            subplot(1,2,modeOFF);
            hold on
            
            plotErrStats(statsOne{modeOFF}{modeON,1},'k',Nexamples);
            plotErrStats(statsOne{modeOFF}{modeOFF,2},'g',Nexamples);
            plotErrStats(statsInf{modeOFF}{modeOFF,2},'b',Nexamples);
            
            % etc.
            title(titlestr{modeOFF});
            axis equal;
            hold off
            
        end
        legend('input (clamped)','one-pass','convergence');
        
        
        
    case 1
        
        % inputvar = statsOnePass{1}.cov + statsOnePass{3}.cov;
        % outputvar = statsFinal{length(params.mods) + 2}.cov;
        
        Nmods = length(params.mods);
        for i = 1:Nmods
            othermods = setdiff(1:Nmods,i);
            
            cov = 0;
            for j = 1:length(othermods)
                cov = cov + statsOne{i}{othermods(j),1}.cov;
            end
            cov
            statsOne{i}{i,2}.cov
            
            cov = 0;
            for j = 1:length(othermods)
                cov = cov + statsInf{i}{othermods(j),1}.cov;
            end
            cov
            statsInf{i}{i,2}.cov
            
            
            pause;
            fprintf('\n');
            
            
        end
        
end


