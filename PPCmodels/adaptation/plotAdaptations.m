function plotAdaptations(SiF,integL0F,DoBF,stepsVWDRAvg,stepsWDRAvg,stepsEMPAvg,params)
% meant to be called from recalibration3.m, it plots the path of the
% adaptations, starting at the original integrated estimates, and migrating
% "backwards" to the unimodal estimates.  This represents the model's
% internal representations moving in the opposite direction.
%
% USAGE:
%   plotAdaptations(xi(:,:,end),IntegL0(:,:,end),thisDoB,...
%       stepsVWDRAvg,stepsWDRAvg,stepsEMPAvg,params)

%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed indexing of integL0F, based on new format for shatL from
%   estStatsCorePP.m; ditto for xiF
% Revised: 09/20/12
%   -functionized, rationalized, bug fixed, etc.
% Created: 09/19/12
%   by JGM
%-------------------------------------------------------------------------%


% init
clc; close all
figure; hold on;
setColors;

% get the first (no adaptation) and final (full adaptation) integ ests
[xInt0,thInt0,xIntF,thIntF] = getFirstAndLastIntegs(integL0F,DoBF,params);

% get the first and last unimodal est---assuming all stims were identical
%   in this batch!!
xx0 = IK2link(SiF(1,:,1),params.roboparams,1)';
th0 = SiF(1,:,2);
xxF = xx0 - mean(xIntF - xInt0);
thF = th0 - mean(thIntF - thInt0);
scatter(xx0(1),xx0(2),[],VIScolor,'o');
scatter(th0(1),th0(2),[],PROPcolor,'o');
scatter(xxF(1),xxF(2),[],VIScolor,'s');
scatter(thF(1),thF(2),[],PROPcolor,'s')

% plot all the integ ests, from start to finish---and the predicted ones!
plotAdptPath(xx0,th0,stepsEMPAvg,[VIScolor;PROPcolor],'-');
plotAdptPath(xx0,th0,stepsVWDRAvg,[VIScolor;PROPcolor],':');
plotAdptPath(xx0,th0,stepsWDRAvg,[VIScolor;PROPcolor],'--');

% some other plots
% meanscatter = @(X,sty)(scatter(mean(X(:,1)),mean(X(:,2)),sty));
% meanscatter(xInt0,'rx');    meanscatter(thInt0,'gx');
% meanscatter(xIntF,'r*');    meanscatter(thIntF,'g*');
% scatter(xInt0(1,:),xInt0(2,:),'r.'); scatter(thInt0(1,:),thInt0(2,:),'g.');
% scatter(xIntF(1,:),xIntF(2,:),'r+'); scatter(thIntF(1,:),thIntF(2,:),'g+');
% plot all the integ ests, from start to finish---and the predicted ones!
% plotAdptPath(xInt0,thInt0,-stepsEMPAvg,'-');
% plotAdptPath(xInt0,thInt0,-80*stepsVWDRAvg,':');
% plotAdptPath(xInt0,thInt0,-0.0155*stepsWDRAvg,'.');
% axis([params.thmin(1) params.thmax(1) params.thmin(2) params.thmax(2)]);

hold off;


% redo in VIS space
figure; hold on;
plotAdptPath2(xx0,th0,stepsEMPAvg,[VIScolor;PROPcolor],'-',params);
plotAdptPath2(xx0,th0,stepsVWDRAvg,[VIScolor;PROPcolor],':',params);
plotAdptPath2(xx0,th0,stepsWDRAvg,[VIScolor;PROPcolor],'--',params);
[xInt0,thInt0,xIntF,thIntF] = getFirstAndLastIntegs2(integL0F,DoBF,params);
xx0 = SiF(1,:,1);
th0 = FK2link(SiF(1,:,2),params.roboparams,1);
xxF = xx0 - mean(xIntF - xInt0);
thF = th0 - mean(thIntF - thInt0);
scatter(xx0(1),xx0(2),[],VIScolor,'o');
scatter(th0(1),th0(2),[],PROPcolor,'o');
scatter(xxF(1),xxF(2),[],VIScolor,'s');
scatter(thF(1),thF(2),[],PROPcolor,'s');
hold off;


% plot trajactories through time
% plotAdptTraj(stepsVWDRAvg,stepsWDRAvg,stepsEMPAvg,params)


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [xInt0,thInt0,xIntF,thIntF] =...
    getFirstAndLastIntegs(integL0F,DoBF,params)

% get modality (logical) indices
iVis = strcmp(params.mods,'Hand-Position');
iProp = strcmp(params.mods,'Joint-Angle');

% init
nCases = size(DoBF,1);

% malloc
xInt0 = zeros(nCases,2);
xIntF = zeros(nCases,2);
thIntF = zeros(nCases,2);

% just read this off
thInt0 = integL0F(:,:,iProp);

% loop through cases and decode or send through the IK (or both) as nec.
for i = 1:nCases
    xInt0(i,:) = IK2link(integL0F(i,:,iVis),params.roboparams,1);
    thing = decoder(DoBF(i,:),params);
    xIntF(i,:) = IK2link(thing(:,1),params.roboparams,1);
    thIntF(i,:) = thing(:,2);
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotAdptPath(x0,th0,steps,colors,typ)

setColors

A = cat(3,[mean(x0,1)' mean(th0,1)'],steps);
cumAdpt = cumsum(A,3);

for iMod = 1:size(cumAdpt,2)
    x = squeeze(cumAdpt(1,iMod,:));
    y = squeeze(cumAdpt(2,iMod,:));
    h = plot(x,y,typ);
    set(h,'Color',colors(iMod,:))
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [xInt0,thInt0,xIntF,thIntF] =...
    getFirstAndLastIntegs2(integL0F,DoBF,params)

% get modality (logical) indices
iVis = strcmp(params.mods,'Hand-Position');
iProp = strcmp(params.mods,'Joint-Angle');

% init
nCases = size(DoBF,1);

% malloc
thInt0 = zeros(nCases,2);
xIntF = zeros(nCases,2);
thIntF = zeros(nCases,2);

% just read this off
xInt0 = integL0F(:,:,iVis);

% loop through cases and decode or send through the IK (or both) as nec.
for i = 1:nCases
    thInt0(i,:) = FK2link(integL0F(i,:,iProp),params.roboparams,1);
    thing = decoder(DoBF(i,:),params);
    xIntF(i,:) = thing(:,1);
    thIntF(i,:) = FK2link(thing(:,2),params.roboparams,1);
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotAdptPath2(x0,th0,steps,colors,typ,params)
% just like plotAdptPath, except that it converts (in a hard-coded way) all
% points into VIS space

A = cat(3,[mean(x0,1)' mean(th0,1)'],steps);
cumAdpt = cumsum(A,3);

for iCase = 1:size(cumAdpt,3)
    cumAdpt(:,1,iCase) = FK2link(cumAdpt(:,1,iCase),params.roboparams,1)';
    cumAdpt(:,2,iCase) = FK2link(cumAdpt(:,2,iCase),params.roboparams,1)';
end

for iMod = 1:size(cumAdpt,2)
    x = squeeze(cumAdpt(1,iMod,:));
    y = squeeze(cumAdpt(2,iMod,:));
    h = plot(x,y,typ);
    set(h,'Color',colors(iMod,:))
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotAdptTraj(stepsVWDRAvg,stepsWDRAvg,stepsEMPAvg,params)

% init
Ndims = params.Ndims;
Nmods = length(params.mods);
nBatches = size(stepsVWDRAvg,3);


% reshape
vwdr = reshape(stepsVWDRAvg,Ndims*Nmods,nBatches)';
wdr = reshape(stepsWDRAvg,Ndims*Nmods,nBatches)';
stepE = reshape(stepsEMPAvg,Ndims*Nmods,nBatches)';

% some useful functions
betas = @(X,Y)(Y'*X)/(X'*X);
aug = @(X)([X,ones(size(X))]);

% loop through data
for i = 1:(Ndims*Nmods)
    figure; hold on;
    
    %%% "empirically determined"---just so you can plot on one graph...
    eta1 = betas(aug(vwdr(:,i)),stepE(:,i));
    eta2 = betas(aug(wdr(:,i)),stepE(:,i));
    
    plot(aug(vwdr(:,i))*eta1','b');
    plot(aug(wdr(:,i))*eta2','r');
    plot(stepE(:,i),'g');
    %     plot(vwdr(:,i)*mean(stepE(:,i))/mean(vwdr(:,i)),'b');
    %     plot(wdr(:,i)*mean(stepE(:,i))/mean(wdr(:,i)),'r');
    %     plot(stepE(:,i),'g');
    
end

hold off;

end




