function testAmatrix(A,S)
% Test to make sure the A matrix is any good
% This is a very weak test, but if your A matrix fails it, it is wrong.


% for a random trial....
iTrial = ceil(158*rand);
Nstates = size(A,2);
Ndims = size(S(iTrial).pos,2);

% get the dynamics...
switch Nstates/Ndims %%% not very elegant
    case 1
        thisX = [S(iTrial).pos];
    case 2
        thisX = [S(iTrial).pos S(iTrial).vel];
    case 3
        thisX = [S(iTrial).pos S(iTrial).vel S(iTrial).acc];
end
thisX = thisX';

% malloc
xhat = zeros(size(thisX));

% simulate forward in time---so you're only checking one-step prediction
for t = 1:(length(S(iTrial).t)-1)
    xhat(:,t+1) = A*thisX(:,t);
end

% compare predicted and actual accelerations
figure(103); clf; hold on;
scatter(xhat(end-1,:),xhat(end,:))
scatter(thisX(end-1,:),thisX(end,:),[],'r.')
hold off;

end