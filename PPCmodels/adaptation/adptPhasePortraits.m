function adptPhasePortraits
% draw phase portraits for the two state-space adaptation rules, WDR and
% VWDR, using the parameters of your two-joint setup.
%
% This provides most of the figures for yr SfN 2012 poster.

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -x0 -> s0, with dimension changes
% Revised: 10/10/12
%   -fixed up figure formatting and exporting
% Revised: 10/09/12
%   -added options of plotting dshatVis/dt and dshatProp/dt (but still in
%   the space of the discrepancy).
% Revised: 10/08/12
%   -incorporated simulation of the RBM's adaptation phase portrait from
%   recalibration4.m
% Revised: 10/02/12
%   cleaned up, functionized, etc.
% Created: 10/01/12
%   by JGM
%-------------------------------------------------------------------------%

%%%%% TO DO
% -gainvec vs. gainmatrix---should be done w/if rather than commments
% -ditto for axis
%%%%%


% init
clc; close all;
setColors;
THING = 'dscrp';

switch THING
    case 'dscrp'
        colors = [1 0 1; 0 1 0; 1 0 0; 0 0 1; 0 0 0];  % 'mgrbk';
        titlestr = '$(\hat s_{\theta} - \hat s_x)$';
        fignames = {'WDRppSTD','VWDRppSTD','CRAppSTD','RBMppSTD'};
    case 'vis'
        colors = VIScolor;
        titlestr = '$\hat s_{\theta},s_x$';
        fignames = {'WDRppSTDalt','VWDRppSTDalt','CRAppSTDalt','RBMppSTDalt'};
    case 'prop'
        colors = PROPcolor;
        titlestr = '$\hat s_{\theta},\hat s_x$';
        fignames = {'WDRppSTDalt','VWDRppSTDalt','CRAppSTDalt','RBMppSTDalt'};
end

MDL = 'STD';
% MDL = 'ALLGAINS';
NOISY = 0;
RESTART = 0;

% generate shifts shifts
M = 15;
lb = [-0.05 -0.05];
ub = [0.05 0.05];
shftvecs = tworowgrid(lb,ub,M);                  % prop - vis!!! (radians) 
%%% say that |prophat - vishat| < 0.05 radians

% load
switch MDL
    case 'STD'
        load ../results/numhidswts/Std050.mat
         gainvecs = [12 12; 12 18; 15 15; 18 12; 18 18];
        % gainvecs = [12 12; 12 18; 15 15; 18 12; 18 18];
        % gainvecs = [15 15];
        dir = 'phaseportraits\Std050\';
        etaVWDR = 200;
        etaWDR = 0.15;
        etaCRA = 0.15;
        mag = 3;
    case 'ALLGAINS'
        load ../results/StdAllGains.mat
        gainvecs = [3 3; 15 3; 9 9; 15 15; 3 15];
        dir = 'phaseportraits\StdAllGains\';
        etaVWDR = 50;
        etaWDR = 0.15;
        etaCRA = 0.15;
        mag = 2;
    otherwise
        error('unrecognized RBM model -- jgm\n');
end


% RBM params
params.Ncases = 100;
params.smpls = 15;
nSubjects = 100;
nBatches = 1;


% loop through gain combos
for jj = 1:size(gainvecs,1)
    
    % RULE-BASED ADAPTATION NEEDS THESE
    if NOISY
        [SigmaVis, SigmaProp, x_p, th_p] =...
            getNoisyAdptParams(gainvecs(jj,:),params);
    else
        [SigmaVis, SigmaProp, x_p, th_p] =...
            getCleanAdptParams(gainvecs(jj,:),params);
    end
    
    % MODEL-BASED ADAPTATION
    params.gmin = gainvecs(jj,:);
    params.gmax = gainvecs(jj,:);
    filename = ['emprPP',num2str(params.gmin,'%02.f')];
    if RESTART
        for kk = 1:size(shftvecs,1)
            tic
            shft = shftvecs(kk,:);
            [IntegL0,DoBF,stVWDR,stWDR,stEMP(:,:,kk)] =....
                recalibrationCorePP(nSubjects,nBatches,shft,wts,params);
            toc
        end
        save(filename,'stEMP','shftvecs','M','params')
    else
       load([dir,filename]);
    end

    
    % VWDR: z1dot = -eta*B*z1, where z1 = vishat - prophat
    B = SigmaVis + SigmaProp;
    switch THING
        case 'dscrp'
            wdrVecFld = @(X)(-etaWDR*X);
            vwdrVecFld = @(X)(-etaVWDR*B*X);
            craVecFld = @(X)(-etaCRA*X);
            rbmVecFld = @(X)(mag*squeeze(diff(stEMP,[],2)));
        case 'vis'
            wdrVecFld = @(X)(etaWDR*inv(inv(SigmaVis)+inv(SigmaProp))*inv(SigmaProp)*X);
            vwdrVecFld = @(X)(etaVWDR*SigmaVis*X);
            craVecFld = @(X)(etaCRA*X);
            vLogInd = strcmp(params.mods,'Hand-Position');
            rbmVecFld = @(X)(mag*squeeze(stEMP(:,vLogInd,:)));
        case 'prop'
            wdrVecFld = @(X)(-etaWDR*inv(inv(SigmaVis)+inv(SigmaProp))*inv(SigmaVis)*X);
            vwdrVecFld = @(X)(-etaVWDR*SigmaProp*X);
            craVecFld = @(X)(-etaCRA*X);
            pLogInd = strcmp(params.mods,'Joint-Angle');
            rbmVecFld = @(X)(mag*squeeze(stEMP(:,pLogInd,:)));
    end       
   
    phaseplot(shftvecs,M,colors(jj,:),wdrVecFld,1);
    phaseplot(shftvecs,M,colors(jj,:),vwdrVecFld,2);
    phaseplot(shftvecs,M,colors(jj,:),craVecFld,3);
    phaseplot(shftvecs,M,colors(jj,:),rbmVecFld,4)
    
end


% mabs and legends
xlbl = '$\hat s_{\theta}^1 - \hat s_x^1$';
ylbl = '$\hat s_{\theta}^2 - \hat s_x^2$';
for i = 1:4
    figure(i)
    % if strcmp(THING,'dscrp'), legend(num2str(gainvecs)); end
    xlabel(xlbl,'Interpreter','latex','fontsize',20);
    ylabel(ylbl,'Interpreter','latex','fontsize',20);
    title(['Dynamics of ',titlestr],'Interpreter','latex','fontsize',20);
    set(gcf, 'PaperPositionMode', 'auto');
    axis equal
    axis([-0.05    0.05   -0.05    0.05]);
    % axis([-0.06    0.06   -0.06    0.06]);
    set(gcf,'Units','centimeters','Position',[14 14 14 14])
    
    exportfig(gcf,[fignames{i},'.eps'],'bounds','tight','Color','cmyk');
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function X = tworowgrid(lb,ub,M)

x1 = linspace(lb(1),ub(1),M);
x2 = linspace(lb(2),ub(2),M);
[X1 X2] = meshgrid(x1,x2);
X = [X1(:)'; X2(:)'];

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function phaseplot(X,M,clr,fxn,fignum)

if size(X,1) > size(X,2); X = X'; end
gridify = @(Y,row)(reshape(Y(row,:),M,M));

Xdot = fxn(X);


figure(fignum);
hold on;
h = quiver(gridify(X,1),gridify(X,2),gridify(Xdot,1),gridify(Xdot,2),0);
set(h,'Color',clr);
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [SigmaVis,SigmaProp,x_p,th_p] = getNoisyAdptParams(gains,params)

%%% first just do it for the small space around the center
nBatches = 500;
params.gmin = gains;
params.gmax = gains;
[R,s0,shft] = generatebiaseddata(nBatches,[],params);


% the initial conditions, in prop space
%%% hard-coded indices
x_p = IK2link(s0(1,:,1),params.roboparams,1)';
th_p = s0(1,:,2)';


% get the two posterior covariances
SINSMerr = covInCalc(R,s0,params);
SigmaVis = SINSMerr{1}.cov;
SigmaProp = SINSMerr{2}.cov;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [SigmaVis, SigmaProp, x_p, th_p] = getCleanAdptParams(gains,params)

% set
shft = [0.05; 0.05];
params.gmin = gains;
params.gmax = gains;
Ndims = params.Ndims;

% get Fisher infos
vInd = strcmp(params.mods,'Hand-Position');
pInd = strcmp(params.mods,'Joint-Angle');
FIv = PPCexpectedFI([params.smin(:,vInd) params.smax(:,vInd)],gains(1),params);
FIp = PPCexpectedFI([params.smin(:,pInd) params.smax(:,pInd)],gains(2),params);

% invert into input "error covariances" (in an asymptotic limit)
SILSCerr{1}.cov = inv(FIv);     SILSCerr{1}.mu = zeros(Ndims,1);
SILSCerr{2}.cov = inv(FIp);     SILSCerr{2}.mu = zeros(Ndims,1);

% get the central th and the shifted one
x_p = scalefxn([0.5 0.5],[0;0],[1;1],params.smin(:,pInd),params.smax(:,pInd));
th_p = x_p + shft;
x = [FK2link(x_p,params.roboparams,1)', th_p];

% convert the error covariances into prop space
[SINSCerrMu,SINSCondCov] = SICE(x,params,SILSCerr);
SigmaVis =  squeeze(SINSCondCov(:,:,:,1));
SigmaProp =  squeeze(SINSCondCov(:,:,:,2));

end
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
% plot a trajectory

%%%
% MM = 3000;
% eta = 5;
% z = zeros(Ndims,MM);
% z(:,1) = [-0.03; +0.04];
% for i = 1:MM-1
%     z(:,i+1) = z(:,i) - eta*B*z(:,i);
% end
% plot(z(1,:),z(2,:),'r')

%-------------------------------------------------------------------------%





%     SigmaInteg = inv(inv(SigmaVis) + inv(SigmaProp));
%     sint0 = SigmaInteg*(inv(SigmaVis)*x_p + inv(SigmaProp)*th_p);
%     wdrVecFld = @(X)(-0.1*(X - repmat(sint0,1,M*M)));
%     phaseplot([-0.05+x_p(1),-0.05+x_p(2)],[0.05+x_p(1),0.05+x_p(2)],...
%         M,colors(kk),wdrVecFld,2)
%     phaseplot([-0.05+th_p(1),-0.05+th_p(2)],[0.05+th_p(1),0.05+th_p(2)],...
%         M,colors(kk),wdrVecFld,3)



% figure(2)
% legend(num2str(gainvecs))
% xlabel('$\hat s_x^1$','Interpreter','latex','fontsize',15);
% ylabel('$\hat s_x^2$','Interpreter','latex','fontsize',15);
% title('WDR Phase Portrait (VIS)')


% figure(3)
% legend(num2str(gainvecs))
% xlabel('$\hat s_{\theta}^1$','Interpreter','latex','fontsize',15);
% ylabel('$\hat s_{\theta}^2$','Interpreter','latex','fontsize',15);
% title('WDR Phase Portrait (PROP)')
