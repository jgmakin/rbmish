function [fracWithinOneStdDev, AvgVrnc] = plotKFStuff2(X,Xhat,CvrnMat)
% This function plots the first and second elements X against time, and
% then overlays the distribution of Xhat, using the mean (Xhat) and
% covariance (CrvnMat).  The plot titles assume that these elements
% correspond to (1D) position and velocity, respectively.  The covariance
% is plotted to the 1-std.-dev line.

%-------------------------------------------------------------------------%
% Created: 05/06/13
%   by JGM
%-------------------------------------------------------------------------%


% final time
T = size(X,2);

% plot x position and predicted x position
figure(12);
clf; 
subplot(2,2,1);
hold on;
plot(1:T,X(1,:),'k');
stddevX = squeeze(sqrt(CvrnMat(1,1,:)));
shadedErrorBar(1:T,Xhat(1,:),stddevX,'-g',1);
title('x position vs. time')
xlabel('t');
ylabel('x');
hold off


% plot y position and predicted y position
subplot(2,2,2);
hold on
plot(1:T,X(2,:),'k');
stddevY = squeeze(sqrt(CvrnMat(2,2,:)));
shadedErrorBar(1:T,Xhat(2,:),stddevY,'-g',1);
title('y position vs. time')
xlabel('t');
ylabel('y');
hold off


% plot the 2D position trajectory
subplot(2,2,3);
hold on
plot(X(1,:),X(2,:),'k');
plot(Xhat(1,:),Xhat(2,:),'g');
%%% axis([-66.8838  126.7498  -91.2477   61.4730]);
% axis([-65 175 -17 123]);
% axis([0 100 -17 123]);
title('Position Trajectories')
xlabel('x')
ylabel('y')
hold off


% plot the 2D velocity trajectory
subplot(2,2,4);
hold on
plot(X(3,:),X(4,:),'k');
plot(Xhat(3,:),Xhat(4,:),'g');
title('Velocity Trajectories')
xlabel('dx/dt')
ylabel('dy/dt')
hold off


% see how bad it is....
vrnc = zeros(size(Xhat));
for tt = 1:size(CvrnMat,3)
    vrnc(:,tt) = sqrt(diag(CvrnMat(:,:,tt)));
end
fracWithinOneStdDev = mean((X<(Xhat+vrnc))&(X>(Xhat-vrnc)),2);

% average it across samples
AvgVrnc = mean(vrnc,2);






end