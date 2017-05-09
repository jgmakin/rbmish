function adaptExpl
% ADAPTEXPLR    Adaptation explorer
%   ADAPTEXPLR simulates adaptation to a "shift" (like prism goggles) to
%   one sensor (i.e., modality), using the WDR (blue) and the VWDR (red)
%   (see Ghahramani 1997).  The sensors are related by a linear xformation,
%   A.
%-------------------------------------------------------------------------%
% Revised: 01/20/12
%   -added the linear xformation
%   -made better plotting functions
%   -functionized
% Created: 11/22/12
%   by JGM
%-------------------------------------------------------------------------%

clc; % close all;

% init
m = 2;
N = 100000; % 8000;
eta = 0.00005; % 0.0006; % 0.002;
setColors;

% do
[covVv,covPp,covVp,covPv,shftp,A] = getParams(m);
[wVp,wPp] = getIntegWts(covVp,covPp);

% calculating shifts in prop space
[deltasp,xoptap,xoptbp] = doAdapt(covVp,covPp,shftp,m,eta,N);
plotAdaptations(wVp,wPp,shftp,deltasp);
plotCalibration(deltasp,covVp,covPp,xoptap,xoptbp,shftp,'prop');
[deltasv,xoptav,xoptbv,shftv] = convertResults(deltasp,xoptap,xoptbp,shftp,A);
plotCalibration(deltasv,covVv,covPv,xoptav,xoptbv,shftv,'vis');

% calculating shifts in vis space
[deltasv,xoptav,xoptbv] = doAdapt(covVv,covPv,shftv,m,eta,N);
plotAdaptations(wVp,wPp,shftp,deltasp);
plotCalibration(deltasv,covVv,covPv,xoptav,xoptbv,shftv,'vis');
[deltasp,xoptap,xoptbp,shftp] = convertResults(deltasv,xoptav,xoptbv,shftv,inv(A));
plotCalibration(deltasp,covVp,covPp,xoptap,xoptbp,shftp,'prop');

% some extra things
% plotDiscrep(eoptV,eoptP)
% deriveEqpt


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [covV,covP,covVp,covPv,shft,A] = getParams(m)
% This function assumes that v = A*p

if m == 1
    covP = 4;
    covV = 2;
    shft = 5;
else
    % make up some params
    rho = 0.3; sig1 = 3; sig2 = 5;
    covP = [sig1^2 sig1*sig2*rho; sig1*sig2*rho sig2^2];
    rho = 0.1; sig1 = 5; sig2 = 2;
    covV = [sig1^2 sig1*sig2*rho; sig1*sig2*rho sig2^2];
    shft = [5;-4];

    % (randomish) xformation b/n spaces
    phi = pi/3;
    A = 1.4*[1.3*cos(phi) -sin(phi); sin(phi) 0.5*cos(phi)];
    
    % convert into prop space
    covVp = inv(A)*covV*inv(A)';
    covPv = A*covP*A';
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [wV,wP] = getIntegWts(covV,covP)
% this fxn assumes that the covs are in the same space

% integration weights
wV = covP/(covP + covV);
wP = covV/(covP + covV);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [deltas,xopta,xoptb] = doAdapt(covV,covP,shft,m,eta,N)
% This function assumes that the covariances, covV and covP, are in the
% same space!  Likewise, shft must be in that space as well.

% malloc
xopta = zeros(m,N);
xoptb = zeros(m,N);
% eoptV = zeros(m,N);
% eoptP = zeros(m,N);
deltas = zeros(N,m,4);

% integration weights
[wV,wP] = getIntegWts(covV,covP);

% loop
deltaVa = 0; deltaPa = 0; deltaVb = 0; deltaPb = 0;
for i = 1:N
    
    % stimulus (uniform prior...)
    x = 2*(rand(m,1)-0.5);
    
    % estimates
    % vhat = sqrtm(A*covVp*A')*randn(m,1) + A*x(:);
    % vhatp = inv(A)*vhat;
    
    vhat = sqrtm(covV)*randn(m,1) + x(:);
    phat = sqrtm(covP)*randn(m,1) + x(:)+shft(:);
    % normrnd(x+shft,sqrtm(covP));
    
    % shift by adaptation [a]
    Va = vhat + deltaVa;
    Pa = phat + deltaPa;
    
    % shift by adaptation [b]
    Vb = vhat + deltaVb;
    Pb = phat + deltaPb;
    
    % compute adaptations
    deltaVa = deltaVa + eta*wP*(Pa - Va);
    deltaPa = deltaPa + eta*wV*(Va - Pa);
    
    deltaVb = deltaVb + eta*covV*(Pb - Vb)/20;          %%% /20 to smooth
    deltaPb = deltaPb + eta*covP*(Vb - Pb)/20;
    
    % compute optimal (it should stay fixed, + noise)
    xopta(:,i) = wV*Va + wP*Pa;
    xoptb(:,i) = wV*Vb + wP*Pb;
%     eoptV(:,i) = Va - (wV*Va + wP*Pa);
%     eoptP(:,i) = Pa - (wV*Va + wP*Pa);
    
    deltas(i,:,1) = deltaVa;
    deltas(i,:,2) = deltaPa;
    deltas(i,:,3) = deltaVb;
    deltas(i,:,4) = deltaPb;
    
    
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotAdaptations(wV,wP,shft,deltas)

plotDeltas(deltas(:,:,[1,3]),wP*shft);
plotDeltas(deltas(:,:,[2,4]),-wV*shft);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotDeltas(deltas,finalvalue)

[N,m,p] = size(deltas);

figure; hold on;
title('visual deltas')
for i = 1:m
    plot(deltas(:,i,1));
    plot(deltas(:,i,2),'r');
    plot((1:N),ones(1,N)*finalvalue(i),'k:')
end
legend('WDR','VWDR','\delta(\infty)');
hold off;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotCalibration(deltas,covV,covP,xopta,xoptb,shft,mod)

figure; hold on;
ttlstr = ['summary plot (no time), ',mod,' space'];
title(ttlstr)

% vis
plot(deltas(:,1,1),deltas(:,2,1),'b');                  % WDR
plot(deltas(:,1,3),deltas(:,2,3),'r');                  % VWDR
error_ellipse(covV,[0,0],'style','c');

% prop
plot(deltas(:,1,2)+shft(1),deltas(:,2,2)+shft(2),'b');  % WDR
plot(deltas(:,1,4)+shft(1),deltas(:,2,4)+shft(2),'r')   % VDWR
error_ellipse(covP,shft,'style','y');

% integ
scatter(mean(xopta(1,:)),mean(xopta(2,:)),'b');         % WDR
scatter(mean(xoptb(1,:)),mean(xoptb(2,:)),'r');         % VDWR

axis equal; legend('WDR','VWDR'); hold off

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [deltasv,xoptav,xoptbv,shftv] =...
    convertResults(deltas,xopta,xoptb,shft,A)


deltasv = zeros(size(deltas));
for i = 1:4
    deltasv(:,:,i) = (A*deltas(:,:,i)')';
end
shftv = A*shft;
xoptav = A*xopta;
xoptbv = A*xoptb;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotDiscrep(eoptV,eoptP)

[m,N] = size(eoptV);

figure; hold on;
title('unimodal/optimal discrepancy');
for i = 1:m
    plot(eoptV','Color',VIScolor);
    plot(eoptP','Color',PROPcolor);
    plot(1:N,zeros(1,N),'k:')
end
legend('vis','prop');
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function deriveEqpt
% Derivation of equilibrium 

syms z WWp dd wwp eps1 eps2 etaa
WW = [-wwp wwp; (1-wwp) -(1-wwp)];


% Learning rules (x = deltaVa, y = deltaPa):
% x[t+1] = x[t] + eta*wP*(Pa - Va)
% y[t+1] = y[t] + eta*wV*(Va - Pa)

% Expanding by def'n of Pa and Va:
% x[t+1] = x[t] + eta*wP*[(phat + y[t-1]) - (vhat + x[t-1])]
% y[t+1] = y[t] + eta*wV*[(vhat + x[t-1]) - (phat + y[t-1])]

% Rearranging:
% x[t+1] = x[t] + eta*wP*y[t-1] - eta*wP*x[t-1] + eta*wP*(phat - vhat)
% y[t+1] = y[t] + eta*wV*x[t-1] - eta*wV*y[t-1] + eta*wV*(vhat - phat)

% Organizing into a vector-matrix equation:
% s[t+1] = s[t] + eta*WW*s[t-1] - eps,
%   where WW = [-wP wP; wV -wV)] and
%   eps = eta*[wP*normal(delta,varP + varV),
%              wV*normal(-delta,varP + varV)].

% Z-tranforming:
% z*S[z] - s[0] = S[z] + eta*W*S[z]/z + eps*z/(z-1)

% Also, s[0] = eps.  So:

% S[z] = inv(z*I - I - eta*WW/z)*(eps*z/(z-1) + eps).
Sz = inv(z*eye(2) - eye(2) - etaa*WW/z)*([eps1; eps2]*z/(z-1) + [eps1; eps2]);

% final value theorem: multiply by (z-1) and evaluate limit as z->inf
p = (z-1)*Sz;

% substitute the expected value of the noise terms
p = subs(p,eps1,etaa*wwp*dd);
p = subs(p,eps2,etaa*(1-wwp)*(-dd));

% simplify (otherwise matlab gives a divide-by-zero error)
p = simple(p);

% "take limit" (substitute z = 1)
subs(p,z,1)

% Hence, learning "stops" when Va is on average equal to Pa; i.e., when
% (deltaVa - deltaPa) = shft.  At this point:
%       deltaVa = wP*shft,
%       deltaPa = -wV*shft.

end
%-------------------------------------------------------------------------%
    