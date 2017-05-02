function [m,b,Rsq] = getrfparams(V0,ldsDATA,params)
% getrfparams   ...
% USAGE:
%   [m,b,Rsq] = getrfparams(V0,ldsDATA)

%-------------------------------------------------------------------------%
% Revised: 07/16/14 (JGM)
%   -functionized, "cleaned up"
%   -incorporated use of scalefxn.m
% Created: 06/??/14
%   by BKD
%-------------------------------------------------------------------------%

% params
N = params.N;
Nhid = size(V0,2);
smin = params.smin(strcmp(params.mods,params.NS));
smax = params.smax(strcmp(params.mods,params.NS));
trajmin = smin;
trajmax = N/(N-1)*(smax - smin) + smin;
thresh = 0.5;
TOPLOT = 0;
upperBound = 0.4; %%% NB: YOU HAVE TO SET THIS MANUALLY!!!  
                  %%% Should probably be made into an argument                 
upperBound = 6;
% upperBound = 12;



% malloc
beta1 = zeros(1,Nhid);
beta0 = zeros(1,Nhid);
Rsq = zeros(1,Nhid);


% positions and (lagged) velocities
%%% hard-coded that Ndims = 1
pos = squeeze(ldsDATA.S(:,1,strcmp(params.mods,params.NS),2:end));
vel = squeeze(ldsDATA.Z(:,2,1:end-1));

% for each hidden unit...
for iUnit = 1:Nhid
    
    % responses of hidden units
    act = squeeze(V0(:,iUnit,2:end));

    
    v = vel(act > thresh);
    x = pos(act > thresh);
     
    if isempty(v)
        fprintf('skipping unit %i: never above thresh\n',iUnit);
        beta0(iUnit) = NaN;
        beta1(iUnit) = NaN;
    else
        xScaled = scalefxn(x,trajmin,trajmax,0,2*pi);
        [beta0(iUnit),beta1(iUnit)] = circle_regress(v,xScaled,[0 upperBound]);
        Rsq(iUnit) = getRsq(xScaled,beta1(iUnit)*v+beta0(iUnit),'RHOSQ');
        
        if TOPLOT
            %%% > ~0.2 for RHOSQ
            titlestr = ['unit ', num2str(iUnit),', R^2 = ',...
                num2str(Rsq(iUnit),'%0.3f')];
            plotSpikes(v,x,beta0(iUnit),beta1(iUnit),trajmin,trajmax,titlestr);
            pause();
        end

    end
end

m = (trajmax-trajmin)*beta1/(2*pi);
b = scalefxn(mod(beta0,2*pi),0,2*pi,trajmin,trajmax);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function plotSpikes(x1,x2,beta0,beta1,xmin,xmax,titlestr)

M = 100;
vPts = linspace(min(x1),max(x1),M)';
xhat = scalefxn(mod(beta1*vPts + beta0,2*pi),0,2*pi,xmin,xmax);

scatterhist(x1,x2)
hold all
plot(vPts,xhat,'g.')
ylim([xmin,xmax])
hold off

title(titlestr);
ylabel('position')
xlabel('velocity')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Rsq = getRsq(y,yhat,METHOD)
% compute coefficient of determination (R^2)

circvar = @(th)(1-sqrt((sum(cos(th)))^2 + (sum(sin(th)))^2)/length(th));
circmean = @(th)(sum(exp(1i*th))/length(th));
circcorr = @(th1,th2)(sin(th1)'*sin(th2)/...
    sqrt((sin(th1)'*sin(th1))*(sin(th2)'*sin(th2))));
        

switch METHOD
    case 'VARIANCEBASED'
        e = y - yhat;
        circvarE = circvar(e);
        circvarT = circvar(y);
        Rsq = 1 - circvarE/circvarT;
    case 'VARPLUSMEANSQ'
        e = y - yhat;
        circvarE = circvar(e);
        MSE = circvarE + angle(circmean(e))^2;
        circvarT = circvar(y);
        Rsq = 1 - MSE/circvarT;
    case 'MSEBASED'
        e = y - yhat;
        MSE = angle(circmean(e.^2));
        circvarT = circvar(y);
        Rsq = 1 - MSE/circvarT;
    case 'RHOSQ'
        rho = circcorr(y,yhat);
        cos(y)'*cos(yhat)/sqrt((cos(y)'*cos(y))*(cos(yhat)'*cos(yhat)));
        Rsq = rho^2;
    otherwise
        error('unrecognized option for computing R^2 -- jgm');
end


end
%-------------------------------------------------------------------------%