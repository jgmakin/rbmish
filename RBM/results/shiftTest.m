% response to discrepant inputs

RESTART = 1;
if RESTART
    shiftLoop;
else
    clear; clc;
    
    % load
    load results/new/eStats1515tight % 1515, 2505, 0525, 1515, 2010, 1020,
    
    % init
    Ndims = params.Ndims;
    Nmods = params.Nmods;
    indices = reshape(1:Nmods*Ndims,Ndims,Nmods);
    switch params.NS
        case 'Hand-Position'
            biasstr = 'visbias';
        case 'Joint-Angle'
            biasstr = 'propbias';
        case 'Gaze-Angle'
            biasstr = 'eyebias';
    end
    NSIND = find(strcmp(params.mods,params.NS));
    indout = indices(:,NSIND);
    nRadii = length(shiftvec);
    nDirections = (size(shiftArray,2) - 1)/nRadii;
end



%-------------------------------------------------------------------------%
% plot all
%%
% close all;
% setColors;

nRadii = length(shiftvec);
i = 2;
% i = 1+nDirections*(2-1)+2;
for iShift = 1:nRadii
    for iDirection = 0:nDirections-1
        
        i = iDirection+nDirections*(iShift-1)+2;    
        dispErrCovs(ErrorStatsArray(i,[1,3:4]),params,40000,2*iShift);
        %%%% HACK: writing in 40000 explicitly here is only "correct" b/c
        %%%% you (think you) used nbatches = 1000 as the argument to
        %%%% DATAGENPP inside test.m.  You always do, but this should still
        %%%% be fixed
        
        % for later!
        shft = shiftArray(:,i);
        biasINP(:,i) = ErrorStatsArray{i,1}{NSIND}.mu;
        biasRBM(:,i) = ErrorStatsArray{i,3}{NSIND}.mu;
        biasOPT(:,i) = ErrorStatsArray{i,4}{NSIND}.mu;
        ecovINP(:,:,i) = ErrorStatsArray{i,1}{NSIND}.cov;
        ecovRBM(:,:,i) = ErrorStatsArray{i,3}{NSIND}.cov;
        ecovOPT(:,:,i) = ErrorStatsArray{i,4}{NSIND}.cov;       
                
        i=i+1;
        
       
    end
    
    figname = ['shifts',num2str(shiftvec(iShift)),'std'];
    % saveas(gcf,['results\figs\',figname],'pdf');
end
bexcess = biasRBM - biasOPT;

for i = 3:2:15
    figure(i);
    close
end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% get some useful quantities
shftmag = sqrt(diag(shiftArray'*shiftArray));
excsmag = sqrt(diag(bexcess'*bexcess));
INPbmag = sqrt(diag(biasINP'*biasINP));
RBMbmag = sqrt(diag(biasRBM'*biasRBM));
OPTbmag = sqrt(diag(biasOPT'*biasOPT));

SHFTMAG = reshape(shftmag(2:end),nDirections,nRadii)';
EXCSMAG = reshape(excsmag(2:end),nDirections,nRadii)';
INPBMAG = reshape(INPbmag(2:end),nDirections,nRadii)';
RBMBMAG = reshape(RBMbmag(2:end),nDirections,nRadii)';
OPTBMAG = reshape(OPTbmag(2:end),nDirections,nRadii)';
%-------------------------------------------------------------------------%



% plot |biases|
%%
figure(); hold on;
A = get(gca,'ColorOrder');
orange = [1 0.5 0];
A = [A; orange];
set(gca,'ColorOrder',A);
shiftmat = repmat(shiftvec,nDirections,1)';
for iDirection = 0:nDirections-1
%     scatter(shftmag(iDirection+2:nDirections:end),...
%             excsmag(iDirection+2:nDirections:end))
    scatter(shiftvec,excsmag(iDirection+2:nDirections:end))
        
    phi = 2*pi*iDirection/nDirections;
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    DRCTCELL{iDirection+1} = num2str(round(180*phi/pi));
end
% plot(SHFTMAG,EXCSMAG);
plot(shiftmat,EXCSMAG);
legend(DRCTCELL{:});
title('Excess-Bias Magnitudes (Marginal Errors)');
xlabel('|input discrepancy|');
ylabel('|output bias - optimal bias|');
hold off
figname = 'excessbiases';
saveas(gcf,['results\figs\',figname],'pdf');

% plot |covariances|
%%
figure(); hold on;
A = get(gca,'ColorOrder');
orange = [1 0.5 0];
A = [A; orange];
set(gca,'ColorOrder',A);
for iDirection = 0:nDirections-1
    for iShift = 1:nRadii
        i = (iShift-1)*nDirections + iDirection + 2;
%         detINP(iShift,iDirection+1) = det(ecovINP(:,:,i));
%         detRBM(iShift,iDirection+1) = det(ecovRBM(:,:,i));
%         detOPT(iShift,iDirection+1) = det(ecovOPT(:,:,i));
        [V D] = eig(ecovOPT(:,:,i)\ecovRBM(:,:,i));
        covdist(iShift,iDirection+1) = norm(log(diag(D)));
        
        [V D] = eig(ecovOPT(:,:,i)\ecovINP(:,:,i));
        cacadist(iShift,iDirection+1) = norm(log(diag(D)));
    end
%     scatter(SHFTMAG(:,iDirection+1),...
%         (detRBM(:,iDirection+1) - detOPT(:,iDirection+1))./detOPT(:,iDirection+1));
%     figure(1); hold on
    % scatter(SHFTMAG(:,iDirection+1),covdist(:,iDirection+1));
    scatter(shiftvec,covdist(:,iDirection+1));
    
%     figure(2); hold on
%     scatter(SHFTMAG(:,iDirection+1),cacadist(:,iDirection+1),'x');
end




% plot(SHFTMAG,(detRBM - detOPT)./detOPT);
% figure(1);
% plot(SHFTMAG,covdist);
plot(shiftmat,covdist);

% figure(2);
% plot(SHFTMAG,cacadist);

hold off;
title('RBM/Optimal Error-Covariance Discrepancy');
xlabel('|input discrepancy|');
% ylabel('(|RBM| - |OPT|)/|OPT|');
ylabel('cov discrepancy')
legend(DRCTCELL)
hold off;
figname = 'excesscovs';
saveas(gcf,['results\figs\',figname],'pdf');


% detINPavg = mean(detINP,2);
% detOPTavg = mean(detOPT,2);
% figure; hold on;
% 
% plot(SHFTMAG,detRBM);
% plot(diag(SHFTMAG),detINPavg,'m:');
% plot(diag(SHFTMAG),detOPTavg,'k:');
% foo = DRCTCELL;
% foo{end+1} = 'input';
% foo{end+1} = 'optimal';
% legend(foo);

%-------------------------------------------------------------------------%






% quiver(shiftArray(1,:),shiftArray(2,:),bexcess(1,:),bexcess(2,:));
% quiver(shiftArray(1,:),shiftArray(2,:),biasOPT(1,:),biasOPT(2,:),'g');
% quiver(shiftArray(1,:),shiftArray(2,:),biasRBM(1,:),biasRBM(2,:),'m');
% title('Error Covariances and Biases for Discrepant Inputs')
% xlabel('radians');
% ylabel('radians');
% axis equal; axis tight; hold off;