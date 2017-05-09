function generalization
% "Local error correction with smoothness constraints"
%
% NB that you might want to change this so that it's called by a script
% that loads the wts and params, and possibly sets the shft, etc. etc.  As
% it is, this function takes no arguments.
%
% Takes about 10 min. to run....

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -changed to force longdata everywhere
% Revised: 12/10/13
%   -changed shatNOPT, shatNPRE, shatNPOST indexing to account for change
%   in estStatsCorePP indexing of its (third and) fourth argument(s).
% Created: 10/02/12
%   by JGM
%-------------------------------------------------------------------------%

% load the standard model
load('results\finalwts\wtsStandard140613.mat','wts','params');

% init
shft = [0.1,0.1];
params.smpls = 15;

% get a subset of the original training data
[localmin,localmax,S0,Q0] = getTrainingSubset(200,0.3,params);

% get new weights by training on these shifted data
wtsNew = wts;
for i = 1:1 % 5
    wtsNew = retrainOnShiftedData(S0,Q0,shft,wtsNew,params);
    fprintf('%i \n',i)
end

% generate new (long) data, on a grid
[R0,S0] = datagenGRID(50,params);

% test on the "normal" data...
[~,~,~,shatNOPT]=estStatsCorePP(S0,params,'CoM',R0);% ...optimal integ
shatNPRE = testTheseData(R0,S0,wts,params);         % ...the original wts
shatNPOST = testTheseData(R0,S0,wtsNew,params);     % ...the new wts

% plot
k = find(strcmp(params.mods,params.NS));
plotPostAdptShifts(S0(:,:,k),shatNOPT(:,:,k),...
    shatNPRE(:,:,k),shatNPOST(:,:,k),localmin,localmax);


5
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [localmin,localmax,S0,Q0] = getTrainingSubset(M,datarng,params)

% temporarily change thmin, thmax to trick getStimuliTiled
pInd = strcmp(params.mods,'Joint-Angle');
cntr = scalefxn([0.5 0.5],[0;0],[1;1],params.smin(:,pInd),params.smax(:,pInd));
localmin = cntr - datarng/2;        localmax = cntr + datarng/2;
params.roboparams.thmin = localmin;	params.roboparams.thmax = localmax;

% get MxM uniform stimuli on a small, central grid
[S0,Q0] = getStimuliTiled(M^2,class(datarng),params);

% shuffle them so that the RBM training doesn't see an order
[yy,ii] = sort(rand(size(S0,1),1));
S0 = S0(ii,:,:);

end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function wtsNew = retrainOnShiftedData(S0,Q0,shft,wts,params)

% init
Ncases = params.Ncases;
Nexamples = size(S0,1);

% offset theta by the input discrepancy
Xshft = S0;
Xshft(:,:,strcmp(params.mods,'Joint-Angle')) =...
    Xshft(:,:,strcmp(params.mods,'Joint-Angle')) + repmat(shft,Nexamples,1);

% generate spike counts (population codes)
Rshft = params.getData(Xshft,Q0);
%%% somewhat hacky...
%%% Rshft = newTunerData(Xshft,params.g*ones(1,params.Ndims),params);
%%% Rshft = shortdata(Ncases,3,Rshft);

% retrain on these data
wtsNew = train(shortdata(Ncases,3,Rshft),wts,params);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [R0,S0] = datagenGRID(M,params)

% get a whole set of training data
[S0,Q0] = getStimuliTiled(M^2,'double',params);

% generate spike counts (population codes)
R0 = params.getData(S0,Q0);
%%% R0 = newTunerData(S0,params.g*ones(1,params.Ndims),params);



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function shatN = testTheseData(R0,S0,wts,params)

% get the model's best integrated estimates
R1 = updownDBN(R0,wts,params,'Nsamples');
[~,~, ~, shatN] = estStatsCorePP(S0,params,'CoM',R1);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotPostAdptShifts(th,thhatOPT,thhatPRE,thhatPOST,localmin,localmax)


figure; hold on;

% % use scatter
% % init
% M = 200;
% 
% % scatter
% scatter(th(:,1),th(:,2),'m.')                       % true
% scatter(thhatOPT(:,1),thhatOPT(:,2),'k.')           % optimal
% scatter(thhatPRE(:,1),thhatPRE(:,2),'r.');          % integs PRE adpt
% scatter(thhatPOST(:,1),thhatPOST(:,2),'g.');        % integs POST adpt


% % use a lattice
% figure; hold on
% Y = thhatPRE; clr = 'r'; 
% YA = reshape(Y(:,1),50,50); YB = reshape(Y(:,2),50,50);
% plot(YA,YB,clr); plot(YA',YB',clr);
% Y = thhatPOST; clr = 'g'; 
% YA = reshape(Y(:,1),50,50); YB = reshape(Y(:,2),50,50);
% plot(YA,YB,clr); plot(YA',YB',clr);
% hold off;


% use quiver
P = sqrt(size(th,1));
gridify = @(Y,row)(reshape(Y(row,:),P,P));
quiver(gridify(th',1),gridify(th',2),...
    gridify(thhatPRE'-thhatOPT',1),gridify(thhatPRE'-thhatOPT',1),0,'b')
quiver(gridify(th',1),gridify(th',2),...
    gridify(thhatPOST'-thhatPRE',1),gridify(thhatPOST'-thhatPRE',1),0,'r')

% labels
% legend('true','optimal','PRE','POST')
legend('PRE-OPT','POST-PRE')
xlabel('$\hat s_{\theta}^1$','Interpreter','latex','fontsize',15);
ylabel('$\hat s_{\theta}^2$','Interpreter','latex','fontsize',15);
title('Generalization of Adaptation')

% draw a box around the points you trained on
plot(linspace(localmin(1),localmax(1),M),localmin(2)*ones(M,1),'k');
plot(linspace(localmin(1),localmax(1),M),localmax(2)*ones(M,1),'k');
plot(localmin(1)*ones(M,1),linspace(localmin(2),localmax(2),M),'k');
plot(localmax(1)*ones(M,1),linspace(localmin(2),localmax(2),M),'k');

hold off

axis equal

end
%-------------------------------------------------------------------------%


