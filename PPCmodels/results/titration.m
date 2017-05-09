% titration
%-------------------------------------------------------------------------%
% slowly reduce the number of output neurons and check the effect on the
% (best) trained network's integration.
%
% uses 1D binomial-binomial network
%-------------------------------------------------------------------------%

clear; clc;
load titration.mat

% variances
xtrct = @(x)(x(end-6:end));
vv = xtrct(v);
pp = xtrct(p);
aa = xtrct(a);
bb = xtrct(b);
ee = xtrct(e);
vecvec = xtrct(vec);

% plot
bar(repmat(vecvec',1,2),[vv;pp]');
hold on
y = [2:0.1:18];
plot(y,y*0 + a(1),'b:')
plot(y,y*0 + b(1),'r:')
plot(y,y*0 + e(1),'g:')
hold off
legend({'$\sigma_{v,out}^2$','$\sigma_{p,out}^2$',...
    '$\sigma_{v,in}^2$','$\sigma_{p,in}^2$','$\sigma_{opt}^2$'},...
    'Interpreter','latex');
xlabel('# hidden units')
ylabel('wtserror covariance')

% biases
figure()
yy = [4 6 8 10 12 16];
B = biases(1:length(yy),:);
bar(repmat(yy',1,4),B)