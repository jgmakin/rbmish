function ManyModelsResults(MODEL,sweptParam)
% ManyModelsResults     Plot results for rEFH, EM1, EM2, for many models
%
% USAGE:
%   ManyModelsResults('1DrEFH','Springs')
%   ManyModelsResults('1DrEFH','Dampers')
%
% You have previously trained a rEFHs and Kalman filters (via EM) on a
% series of dynamical model, created by varying one of the three parameters
% of the harmonic oscillator (m, b, k).  These are saved in .mat files.  To
% tell this file which to load, provide a string with the parameters type
% (see above).
 

%-------------------------------------------------------------------------%
% Created: 06/09/15
%   by JGM
%-------------------------------------------------------------------------%



[Nsprings,AllparamsEM1,AllparamsEM2,Allwts,params] =...
    getSpringResults(MODEL,sweptParam);
[eStats,paramsVec] = getSpringStats(Nsprings,sweptParam,...
    AllparamsEM1,AllparamsEM2,Allwts,params);
plotSpringStats(paramsVec,eStats,sweptParam);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Nsprings,AllparamsEM1,AllparamsEM2,Allwts,EFHparams] =...
    getSpringResults(MODEL,sweptParam)


yrdir = 'dynamical\finalwts\';

prefix = 'wts';
load([yrdir,prefix,MODEL,'Many',sweptParam,'.mat'])
EFHparams = params;

prefix = 'LDSparamsEM1stOrd';
load([yrdir,prefix,MODEL,'Many',sweptParam,'.mat'])
AllparamsEM1 = Allparams;

prefix = 'LDSparamsEM2ndOrd';
load([yrdir,prefix,MODEL,'Many',sweptParam,'.mat'])
AllparamsEM2 = Allparams;

Nsprings = length(AllparamsEM1);
if length(AllparamsEM2)~=Nsprings, error('wtf? -- jgm'); end
if length(Allwts)~=Nsprings, error('wtf? -- jgm'); end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [eStats,paramsVec] = getSpringStats(Nsprings,sweptParam,...
    AllparamsEM1,AllparamsEM2,Allwts,params)

% The constants for the A matrix
switch sweptParam
    case 'Springs'
        bs = 0.25*ones(1,Nsprings);
        ks = linspace(0,5,Nsprings);
        ms = 5*ones(1,Nsprings);
        paramsVec = ks;
    case 'Dampers'
        bs = linspace(0,1,Nsprings);
        ks = 3*ones(1,Nsprings);
        ms = 5*ones(1,Nsprings);
        paramsVec = bs;
    case 'Masses'
        bs = 0.25*ones(1,Nsprings);
        ks = 3*ones(1,Nsprings);
        ms = linspace(1,10,Nsprings);
        paramsVec = ms;
end
dt=0.05;

% loop: each "experiment" each of which uses a different spring constant
for iSpring = 1:Nsprings
    params.dynamics.A = [1.0000, dt;...
        -ks(iSpring)/ms(iSpring)*dt, -(bs(iSpring)/ms(iSpring)*dt-1)];
    LDSdataTest = getLDSdata(params);
    
    name = ['EM$^',num2str(size(AllparamsEM1(iSpring).A,2)),'$'];
    pEM = KF4PPC(LDSdataTest,AllparamsEM1(iSpring),name);
    eStats(iSpring,1) = testDynamics(LDSdataTest,params,0,pEM);
    
    name = ['EM$^',num2str(size(AllparamsEM2(iSpring).A,2)),'$'];
    pEM = KF4PPC(LDSdataTest,AllparamsEM2(iSpring),name);
    eStats(iSpring,2) = testDynamics(LDSdataTest,params,0,pEM);
    
    [~,pSENSORY,pEFH] = EFHfilter(LDSdataTest,Allwts{iSpring},params);
    eStats(iSpring,3) = testDynamics(LDSdataTest,params,0,pEFH);
    eStats(iSpring,4) = testDynamics(LDSdataTest,params,0,pSENSORY);
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotSpringStats(paramsVec,eStats,sweptParam)

% plot
figure(1546); clf; hold on;
legendcell = {};
for iModel = 1:(size(eStats,2)-1)
    plot(paramsVec,[eStats(:,iModel).Cvrn] + [eStats(:,iModel).Xpct].^2);
    legendcell = {legendcell{:},eStats(1,iModel).tags.name};
end
axis([paramsVec(1),paramsVec(end),0,max([eStats(:,1).Cvrn])]);
% legend(legendcell)
% xlabel('spring constant (kg/s^2)')
ylabel('MSE')
hold off;

matlab2tikzWrapper(['MSEsMany',sweptParam],figure(1546))

end
%-------------------------------------------------------------------------%