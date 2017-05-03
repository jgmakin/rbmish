function plotKFStuff(X,Xhat,CvrnMat)
% This function plots the first and second elements X against time, and
% then overlays the distribution of Xhat, using the mean (Xhat) and
% covariance (CrvnMat).  The plot titles assume that these elements
% correspond to (1D) position and velocity, respectively.  The covariance
% is plotted to the 1-std.-dev line.

%-------------------------------------------------------------------------%
% Created: 04/30/13
%   by JGM
%-------------------------------------------------------------------------%

% final time
T = size(X,2);

% plot position and predicted position
figure(10);
clf
subplot(1,2,1);
hold on
plot(1:T,X(1,:),'k');
stddevP = squeeze(sqrt(CvrnMat(1,1,:)));
shadedErrorBar(1:T,Xhat(1,:),stddevP,'-g',1);
title('Position vs. Time')
hold off


% plot velocity and predicted velocity
subplot(1,2,2);
hold on
plot(1:T,X(2,:),'k');
stddevV = squeeze(sqrt(CvrnMat(2,2,:)));
shadedErrorBar(1:T,Xhat(2,:),stddevV,'-g',1);
title('Velocity vs. Time')
hold off


% make sure that about 68% of the data lie within the first stddev (this
% requires a long T, ~200)
mean((X<(Xhat+[stddevP';stddevV']))&(X>(Xhat-[stddevP';stddevV'])),2)



end