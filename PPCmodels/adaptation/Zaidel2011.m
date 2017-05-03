% a little reproduction of the simulation in Zaidel2011, Figure 9.
%-------------------------------------------------------------------------%
% Created: 01/06/12
%   by JGM
%-------------------------------------------------------------------------%


clear; clc; close all;

N = 2000;

mu1 = [1; -1];
mu2 = [-1; 1];


for i = 3:-1:1
    
    % new variance
    sig = [i*0.5; 1];
    
    % sample the dstrbs
    r1 = repmat(mu1,1,N) + repmat(sig,1,N).*randn(2,N);
    r2 = repmat(mu2,1,N) + repmat(sig,1,N).*randn(2,N);
    
    % type II regression
    x = [r1(1,:) r2(1,:)];
    y = [r1(2,:),r2(2,:)];
    [m,b,r,sm,sb,xbar,ybar] = lsqfitma(x,y);
    
    % plot samples
    figure; hold on
    scatter(r1(1,:),r1(2,:));
    scatter(r2(1,:),r2(2,:));
    
    % plot regression lines
    xx = min(x):0.1:max(x);
    yy = m*xx + b;
    plot(xx,yy);
    axis equal; hold off;
    
    
end

