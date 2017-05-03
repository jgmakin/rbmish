function cmprAdptTrajs
% for comparing the rates of *prop* adaptation from different runs of
% recalibration3.m, presumably each with different settings of the gains
% (and hence posterior covariances).  If the model adapts like the WDR or
% VWDR, prop should adapt faster for larger prop covariance, i.e. smaller
% prop gain.  However, the WDR predicts that gains = [12 12] will look like
% gains = [15 15] (relative gains matters), whereas the VWDR predicts that
% it will look like gains = [18 12] (only prop gain matters).
%
% The "trajectories" plotted here are just the adaptations steps taken at 
% each point "in time."
%
% NB: Eventually, the noise pushes the exponential approach above zero, in
% which case log(-x) is complex.  J = 40 cuts off the curve before that pt.

%-------------------------------------------------------------------------%
% Created: 09/27/12
%   by JGM
%-------------------------------------------------------------------------%

% init
close all; clc
propind = 2;
colors = 'mgrbk';
% gainvecs = [12 12; 18 12; 15 15; 18 18; 12 18];
gainvecs = [03 03; 03 15; 09 09; 15 03; 15 15];
ext = 'slowWIDE'; % 'slow'; % 'quick'; % 'avg';
% gainvecs = [03 15; 09 09; 15 15];
% gainvecs = [3 15; 3 15; 9 9; 15 3; 15 15];

% load and plot
for i = 1:size(gainvecs,1)
    filename = ['adaptationPaths/adpt',num2str(gainvecs(i,1),'%02.f'),...
        num2str(gainvecs(i,2),'%02.f'),ext];
    load(filename);
    

    figure(1); 
    trajs = squeeze(stEMP(:,propind,:));
    J = 30;
    plotTraj(trajs,colors(i))
    plotWithRgrn(J,trajs,colors(i))
    
    figure(2);
    trajs = squeeze(stVWDR(:,propind,:));
    J = 600;
    plotTraj(trajs,colors(i))
    plotWithRgrn(J,trajs,colors(i))
    
    figure(3);
    J = 600;
    trajs = squeeze(stWDR(:,propind,:));
    plotTraj(trajs,colors(i))
    plotWithRgrn(J,trajs,colors(i))
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function plotTraj(trajs,clr)

K = size(trajs,2);
subplot(2,2,1); 
hold on;
plot(1:K,trajs(1,:),clr)
hold off;

subplot(2,2,2);
hold on;
plot(1:K,trajs(2,:),clr)
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotWithRgrn(J,trajs,clr)

subplot(2,2,3);
plotExpRgrn(trajs(1,1:J)',clr);
    
subplot(2,2,4); 
plotExpRgrn(trajs(2,1:J)',clr);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotExpRgrn(vec,clr)

J = length(vec);
beta = linregress([1:J;ones(1,J)]',log(-vec));
hold on
plot(1:J,-exp(beta(2))*exp(beta(1)*(1:J)),[clr,':']);
plot(1:J,vec,clr);
hold off

end
%-------------------------------------------------------------------------%