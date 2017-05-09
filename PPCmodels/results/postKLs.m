% postKLs
%   Kind of a misnomer.  This script computes and plots (in 3D) the
%   fractional information loss suggested by JMB.

%-------------------------------------------------------------------------%
% Revised: 07/08/14
%   -altered to accommodate new getSuffStats.m
%   -eliminated uses of tex.m
%   -found a bug!!  You had left out a factor of 1/2 in your calculation of
%   the entropy of a Gaussian.  Restoring it increases the maximal
%   information loss---from 1.15% to about 2%.
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

%%%%% TO DO %%%%%%
% (1) parallelize getKLdivs
% (2) parallelize nrmlentropy?
% (3) tikzify figure
%%%%%%%%%%%%%%%%%%


% p = conditioned on inputs (optimal)
% q = conditioned on hidden samples
% r = conditioned on hidden means
% s = conditioned on mean input spike counts, etabar_v and etabar_p


clear; clc;
% load results/finalwts/Std050.mat
load('results\finalwts\wtsStandard140613','wts','params');


% init
setColors;
params.smpls = 15;
M = 10;
nrmlEntropy = @(pcn)(log((2*pi*exp(1))^size(pcn,1)/det(pcn))/2);
unifEntropy = @(smin,smax)(log(prod(smax-smin)));


% init params
Ndims = params.Ndims;
thmin = params.roboparams.thmin;
thmax = params.roboparams.thmax;
Nexamples = 4000;


%%%%%%%%%%%%%
% decode better eta(v)
% N = M^2;
% trainvec = 1:N/2; testvec = trainvec + N/2;
% [etaHATq Rsqq] = improveEta(pVIS,pPROP,qVIS,qPROP,trainvec,testvec); % smpls
% [etaHATr Rsqr] = improveEta(pVIS,pPROP,rVIS,rPROP,trainvec,testvec); % means
%%%%%%%%%%%%%


% malloc
avgKL1 = zeros(M,M,3);                      % with means
avgKL2 = zeros(M,M,3);                      % "without means"

% loop through gains
tic
gains = linspace(mean(params.gmin),mean(params.gmax),M);
Hprop = unifEntropy(thmin,thmax);
for gvind = 1:M
    for gpind = 1:M
        
        % generate data
        params.gmin = [gains(gvind) gains(gpind)];
        params.gmax = [gains(gvind) gains(gpind)];
        R = generateData(Nexamples,params);
        
        % get the different posterior distributions
        [p,q,r] = getEFHposteriors(R,wts,params);
        pINT.Xpct = p.Xpct(:,:,strcmp(p.srcs,'opt'));
        pINT.Info = p.Info(:,:,:,strcmp(p.srcs,'opt'));
        qINT.Xpct = q.Xpct(:,:,strcmp(p.srcs,'opt'));
        qINT.Info = q.Info(:,:,:,strcmp(p.srcs,'opt'));
        rINT.Xpct = r.Xpct(:,:,strcmp(p.srcs,'opt'));
        rINT.Info = r.Info(:,:,:,strcmp(p.srcs,'opt'));
        
        % get KL divergences
        KL1 = getKLdivs(pINT,0,qINT,rINT);
        KL2 = getKLdivs(pINT,1,qINT,rINT);
        %%% interesting that you use qINT rather than qPROP
        
        % get (approx) "baseline" KL divergence, KL{p(s|r)||p(s)}
        condMI = zeros(length(pINT),1);
        for iExample = 1:size(pINT.Info,1)
            condMI(iExample) = Hprop - nrmlEntropy(...
                squeeze(pINT.Info(iExample,:,:)));
        end
        KL1(:,3) = condMI;
        KL2(:,3) = condMI;
        
        
        % store the averages over all trials
        avgKL1(gvind,gpind,:) = mean(KL1,1);
        avgKL2(gvind,gpind,:) = mean(KL2,1);
        
    end
end
toc


% print
fontsize = 25;
[gv,gp] = meshgrid(gains,gains);
h = figure; surf(gv,gp,avgKL1(:,:,1)./avgKL1(:,:,5));
% xlabel('$\gain_{\vis}$','Interpreter','Latex','Fontsize',fontsize)
% ylabel('$\gain_{\prop}$','Interpreter','Latex','Fontsize',fontsize)
% zlabel('normalized KL divergence')
% saveas(h,'C:\#code\RBM\results\postcov\normalizedKLvsGains','eps');
%%% tikzify me




    
% figure;
% grandAvgKL1bits = squeeze(mean(mean(avgKL1,2),1))/log(2);
% grandAvgKL2bits = squeeze(mean(mean(avgKL2,2),1))/log(2);
% h = barstar([grandAvgKL2bits(3),grandAvgKL2bits(1)],[],'bits',...
%     {'p(\stim|\bar\eta^{\vis},\bar\eta^{\prop})','\RBMposterior'},[],...
%     [OPTcolor; NNcolor],'KL divergences from $\OPTposteriorExpl$',...
%     fontsize,[]);
% saveas(h,'C:\#code\RBM\results\postcov\postKLsNoMeansBits','eps');
