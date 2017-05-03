function cmprAdptPaths
% for comparing the paths of (vis-prop) adaptation from different runs of
% recalibration3.m, each with different gains (and hence posterior covs).
% To see that this is the right thing to plot (at least for the WDR and
% VWDR), see your TeX notes.

%-------------------------------------------------------------------------%
% Revised: 10/03/12
%   -changed to plot evolution of diffs rather than evolution of prop
%   steps.
% Created: 10/02/12
%   by JGM
%-------------------------------------------------------------------------%

% init
close all; clc
colors = 'mgrbk';
ext = 'slowWIDE'; % 'slow'; % 'quick'; % 'avg';

gainvecs = [03 03; 15 03; 03 15; 09 09; 15 15];
% gainvecs = [12 12; 18 12; 15 15; 18 18; 12 18];
% gainvecs = [03 15; 09 09; 15 15];
% gainvecs = [3 15; 3 15; 9 9; 15 3; 15 15];

%%%%%
% x = FK2link(scalefxn([0.5;0.5],[0;0],[1;1],[-pi/2; pi/4],[pi/4; 3*pi/4]),params,1);
shft = [0.05; 0.05];
% th = scalefxn([0.5;0.5],[0;0],[1;1],[-pi/2; pi/4],[pi/4; 3*pi/4]) + shft;
getTraj = @(T)(squeeze(diff(cumsum(T,3),[],2)) + repmat(shft,1,size(T,3)));
%%%%%

% load and plot
for i = 1:size(gainvecs,1)
    filename = ['adaptationPaths/adpt',num2str(gainvecs(i,1),'%02.f'),...
        num2str(gainvecs(i,2),'%02.f'),ext];
    load(filename);
   
    %%% these are plots of the difference b/n shatVis and shatProp
    figure(1); 
    trajs = getTraj(stEMP);
    plotWithRgrn(shft,trajs,colors(i))
    plotAdptPath(trajs,colors(i));   
    
    figure(2);
    trajs = getTraj(stVWDR);
    plotWithRgrn(shft,trajs,colors(i))
    plotAdptPath(trajs,colors(i));   
    
    figure(3);
    trajs = getTraj(stWDR);
    plotWithRgrn(shft,trajs,colors(i))
    plotAdptPath(trajs,colors(i));       
    
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
function plotWithRgrn(shft,trajs,clr)

subplot(2,2,1);
plotExpRgrn(trajs(1,:)',shft(1),clr);
    
subplot(2,2,2); 
plotExpRgrn(trajs(2,:)',shft(2),clr);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotExpRgrn(z,z0,clr)

J = length(z);
beta = linregress([1:J]',log(z0./z));
hold on
plot(1:J,z0*exp(-beta(1)*(1:J)),[clr,':']);
plot(1:J,z,clr);
hold off

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotAdptPath(trajs,clr)

subplot(2,2,3);
hold on;  
plot(trajs(1,:),trajs(2,:),clr); 
axis equal; 
hold off;    
    
end
%-------------------------------------------------------------------------%


















