function KFmaster
%   Creates some data for a filtering problem, then runs the KF on it, and
%   plots the results.

%%% To do: add controls!!

%-------------------------------------------------------------------------%
% Revised: 04/30/13
%   -functionized
%   -made it work
% Created: 04/29/13
%   by JGM
%-------------------------------------------------------------------------%

clear; clc; % close all


params = initModel;
[X,Y] = createData(params);
KFdstrbs = KalmanFilter(params,Y);
plotKFStuff(X,KFdstrbs.XHATMU,KFdstrbs.CVRNMU);


checkSteadyState(KFdstrbs.CVRNMU,params)

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function params = initModel

% params
params.T = 100;
dT = 1;

% prior distribution
params.mu0 = [0; 0];
params.Info0 = [0 0; 0 0.3];

% state equations
params.A = [1 dT; 0 1];
params.C = [1 0];


% SigmaX = 5;
params.SigmaX = [4 0; 0 0];
params.SigmaY = 2.5;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [X,Y] = createData(params)

% params
T = params.T;
A = params.A;
C = params.C;

mu0 = params.mu0;
Info0 = params.Info0;

Nstates = size(A,1);
Nobsvs = size(C,1);

if det(params.SigmaX)==0, SigmaXonehalf = sqrt(params.SigmaX);
    %%% this is rather hacky....
else SigmaXonehalf = sqrtm(params.SigmaX); 
end

if det(params.SigmaY)==0, SigmaYonehalf = sqrt(params.SigmaY);
    %%% this is rather hacky....
else SigmaYonehalf = sqrtm(params.SigmaY); 
end


% malloc
X = zeros(Nstates,T);
Y = zeros(Nobsvs,T);

% initialize
if det(Info0)~=0
    X(:,1) = mu0 + sqrtm(inv(Info0))*randn(Nstates,1);
else
    %%% a hack, premised on the zeroed info being on position
    X(1,1) = mu0(1);
    X(2,1) = mu0(2) + randn/sqrt(Info0(2,2));
end
 
% loop
for t = 1:(T-1)
    % x(:,t+1) = A*x(:,t) + sqrtm(SigmaX)*randn(Nstates,1);
    X(:,t+1) = A*X(:,t) + SigmaXonehalf*randn(size(SigmaXonehalf,1),1);
    Y(:,t) = C*X(:,t) + SigmaYonehalf*randn(Nobsvs,1);
end
Y(:,T) = C*X(:,T) + SigmaYonehalf*randn(Nobsvs,1);
fprintf('\n\nData created...\n\n\n')

end
%-------------------------------------------------------------------------%


