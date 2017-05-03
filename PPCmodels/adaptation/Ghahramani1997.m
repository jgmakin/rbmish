function Ghahramani1997
%%%%% simulate ZG's results in Ghahramani1997

%-------------------------------------------------------------------------%
% Revised: 10/07/12
%   -bug fixes
% Created: 10/05/12
%   by JGM
%-------------------------------------------------------------------------%

clear; clc; close all;
setColors

% 1D arc of space, in degrees
M = 100;
x = linspace(-45,45,M);
var = 2;

% variances (constructed to look like Ghahramani1997)
VarVis = 0.05*((x/2).^2) + 0.5;
VarAud = 0.05*( x   .^2) + 10;

% optimal integration weights
WVis = (1./(1./VarVis)+(1./VarAud)).*(1./VarVis);
WAud = (1./(1./VarVis)+(1./VarAud)).*(1./VarAud);


% plot variance (cf. with Fig. 2 in Ghahramani1997)
% figure; hold on;
% plot(x,VarVis,'g');
% plot(x,VarAud,'r');
% hold off


% function to get the optimal estimate
optcombo = @(vis,aud)(((VarVis.^-1 + VarAud.^-1).^-1).*...
    (1./VarVis.*vis + 1./VarAud.*aud));


% loop
delta0 = 15;
N = 1000;
VIS0 = x;       
AUD0 = VIS0 - delta0;
shatOPT = optcombo(VIS0,AUD0);


% WDR
etaWDR = 0.01; % 0.04;
stepfxn = @(delta,deg)(-etaWDR*[WAud(deg)*delta(:,deg); -WVis(deg)*delta(:,deg)]);
VIS = VIS0;     AUD = AUD0;
for i = 1:N
    del = VIS - AUD;
    dgr = ceil(M*rand);
    step = stepfxn(del,dgr);
    VIS = VIS + rbf(step(1),x,dgr,var);
    AUD = AUD + rbf(step(2),x,dgr,var);
end

figure; hold on;
plot(x,VIS-VIS0,'g');
plot(x,AUD-AUD0,'r');
% plot(x,VIS-AUD,'b');
plot(x,shatOPT-VIS0,'g:');
plot(x,shatOPT-AUD0,'r:');
title('WDR Adaptation')
legend('VIS','AUD','final VIS','final AUD')
hold off;


% VWDR
etaVWDR = 0.001; % 0.001;
stepfxn = @(delta,deg)(-etaVWDR*[VarVis(deg)*delta(:,deg); -VarAud(deg)*delta(:,deg)]);
VIS = VIS0;     AUD = AUD0;
caca = 0.1;
for i = 1:N
    del = VIS - AUD;
    dgr = ceil(M*rand);
    step = stepfxn(del,dgr);
    VIS = VIS + rbf(step(1),x,dgr,var);
    AUD = AUD + rbf(step(2),x,dgr,var);
end

figure; hold on;
plot(x,VIS-VIS0,'g');
plot(x,AUD-AUD0,'r');
% plot(x,VIS-AUD,'b');
plot(x,shatOPT-VIS0,'g:');
plot(x,shatOPT-AUD0,'r:');
title('VWDR Adaptation');
legend('VIS','AUD','final VIS','final AUD')
hold off;



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function basisfxn = rbf(wt,x,mu_ind,var)

basisfxn = wt*var*sqrt(2*pi)*normpdf(x,x(mu_ind),var);

end
%-------------------------------------------------------------------------%





