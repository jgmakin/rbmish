function KL = posteriorRecovery(R,wts,params)
% posteriorRecovery    mean squared errors for the gain proxies
%   GAINERRORS computes and plots (&c.) the R^2s for the gain
%   proxies---etas for the input layer, some other thing for the hidden
%   layer.
%
% USAGE:
%   load('results\finalwts\wtsStandard140613.mat')
%   [R0,S] = generateData(1000,params);
%   KL = posteriorRecovery(R0,wts,params);

%-------------------------------------------------------------------------%
% Revised: 10/22/14
%   -now passes strings of color names to tikzBarGraph
% Revised: 07/07/14
%   -major changes, centered around accommodating the new getSuffStats.m
% Revised: 10/24/12
%   -massively
% Adapted: 10/23/12
%   -from gainerrors.m, with huge changes
%   by JGM
%-------------------------------------------------------------------------%


% init
setColors;
params.smpls = 15;
[Ncases,Nvis,Nbatches] = size(R);
Nexamples = Ncases*Nbatches;
trainvec = 1:Nexamples/2; testvec = trainvec + Nexamples/2;

% get the different posterior distributions
[p,q,r] = getEFHposteriors(longdata(R),wts,params);
[pVIS,pPROP,pINT] = separatePosteriors(p);
[qVIS,qPROP,qINT] = separatePosteriors(q);
[rVIS,rPROP,rINT] = separatePosteriors(r);
sINT.Xpct = p.Xpct(:,:,strcmp(p.srcs,'optetabar'));
sINT.Info = p.Info(:,:,:,strcmp(p.srcs,'optetabar'));


% get KL divergences
KL1 = getKLdivs(pINT,0,pVIS,pPROP,qINT,rINT,sINT);  % don't ignore means
KL2 = getKLdivs(pINT,1,pVIS,pPROP,qINT,rINT,sINT);  % ignore means!


% plot KL divs for both single-pop. posteriors, the RBM post., and OPT
avgKLs = [mean(KL1(:,1:end-1),1),mean(KL2(:,1:end-1),1)];
xaxislabels = {'$p(\stims|{\obsvs}^{\visl})$',...
    '$p(\stims|{\obsvs}^{\prop})$','$\RBMposterior$',...
    '$q(\stims|\bar{\ltnts})$','$p(\stim|{\obsvs}^{\visl})$',...
    '$p(\stims|{\obsvs}^{\prop})$','$\RBMposterior$',...
    '$q(\stims|\bar{\ltnts})$'};
titlestr = 'KL divergences from $\OPTposteriorExpl$';
avgKLerrorBars = zeros(length(avgKLs),2);
tikzBarGraph(1:size(avgKLs,2),avgKLs',avgKLerrorBars,0,...
    xaxislabels',[],'nats',titlestr,...
    {'vclr','pclr','rbmclr','rbmclr','vclr','pclr','rbmclr','rbmclr'}',...
    2,0.32,'overlapping',{},...
    ['posteriorKLthing',date]);


% plot KL divs for the RBM posts. and one that gets access to etabar(1,2)
% avgKLsbits = [mean(KL1(:,3:end),1),mean(KL2(:,3:end),1)]/log(2);
avgKLsbits = mean(KL2(:,3:end),1)/log(2);
xaxislabels = {'$q(\stims|\ltnts,\hat\stims(\obsvs))$',...
    '$q(\stims|\bar{\ltnts},\hat\stims(\obsvs))$',...
    '$p(\stims|\bar\ttlspks^{\visl},\bar\ttlspks^{\prop},\hat\stims(\obsvs))$'};
titlestr = 'KL divergences from $\OPTposteriorExpl$';
avgKLerrorBars = zeros(length(avgKLsbits),2);
tikzBarGraph(1:size(avgKLsbits,2),avgKLsbits',avgKLerrorBars,0,...
    xaxislabels',[],'bits',titlestr,...
    {'rbmclr','rbmclr','optclr'}',...
    2,0.62,'overlapping',{},...
    ['posteriorKLdivsbits',date]);


% decode better eta(v)
[ttlSpksHatq,Rsqq] = improveEta(p,q,trainvec,testvec); % smpls
[ttlSpksHatr,Rsqr] = improveEta(p,r,trainvec,testvec); % means
xaxislabels = {'$\ttlspks^{\visl}(\ltnts)$','$\ttlspks^{\prop}(\ltnts)$',...
    '$\ttlspks^{\visl}(\bar\ltnts)$','$\ttlspks^{\prop}(\bar\ltnts)$'};
titlestr = 'Recovery of total spike counts';
RSQerrorBars = zeros(length([Rsqq Rsqr]),2);
tikzBarGraph(1:(size(Rsqq,2)+size(Rsqr,2)),[Rsqq Rsqr]',RSQerrorBars,0,...
    xaxislabels',[],'$R^2$',titlestr,...
    {'rbmclr','rbmclr','rbmclr','rbmclr'}',...
    2,0.32,'overlapping',{},...
    ['etarecovery',date]);


%%% Can even better results be achieved by finding better etas(R1)?
% update the variances based on these improved etas
qINT.Info = updatePrcns(q,ttlSpksHatq,testvec);
rINT.Info = updatePrcns(r,ttlSpksHatr,testvec);

% get KL divergences
pINT.Info = pINT.Info(testvec,:,:); pINT.Xpct = pINT.Xpct(testvec,:);
pVIS.Info = pVIS.Info(testvec,:,:); pVIS.Xpct = pVIS.Xpct(testvec,:);
pPROP.Info = pPROP.Info(testvec,:,:); pPROP.Xpct = pPROP.Xpct(testvec,:);
qINT.Xpct = qINT.Xpct(testvec,:); 
rINT.Xpct = rINT.Xpct(testvec,:);
sINT.Info = sINT.Info(testvec,:,:); sINT.Xpct = sINT.Xpct(testvec,:);

KL3 = getKLdivs(pINT,0,pVIS,pPROP,qINT,rINT,sINT);
mean(KL3,1)

KL4 = getKLdivs(pINT,1,pVIS,pPROP,qINT,rINT,sINT);
mean(KL4,1)

KL = [KL1; KL2; KL3; KL4];


% plot some example posterior covariances
%%% plotPostCovs(pINT,pVIS,pPROP,qINT,rINT);
%%% you broke this (07/07/14, jgm)

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [YHAT,b,MSE,Rsq] = linearfits2(Y,X,trainvec,testvec)
%%% think about returning *all* the YHAT, even though you only do the Rsq
%%% calculation on the held-out data
%%% Also: why not do LOOCV???


% init
M = length(trainvec);
N = length(testvec);

% fit an affine ("linear") model
beta1 = linregress([X(trainvec,1) ones(M,1)],Y(trainvec,1));
beta2 = linregress([X(trainvec,2) ones(M,1)],Y(trainvec,2));
YHAT = [[X(testvec,1) ones(N,1)]*beta1,[X(testvec,2) ones(N,1)]*beta2];
b = [beta1 beta2];

% get error stats
E = Y(testvec,:) - YHAT;
MSE = mean(E.^2);
MST = var(Y(testvec,:));
Rsq = 1 - MSE./MST;


end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [YHAT,b,MSE,RsqCV] = linearfits3(Y,X)
%%% think about returning *all* the YHAT, even though you only do the Rsq
%%% calculation on the held-out data
%%% Also: why not do LOOCV???


% init
N = size(X,1);
X1 = [X(:,1) ones(N,1)];        Y1 = Y(:,1);
X2 = [X(:,2) ones(N,1)];        Y2 = Y(:,2);

% fit an affine ("linear") model
[beta1, RsqCV1] = linregress(X1,Y1,'LOO');
[beta2, RsqCV2] = linregress(X2,Y2,'LOO');

% get error stats
YHAT = [X1*beta1,X2*beta2];
E = Y - YHAT;
MSE = mean(E.^2);

% collect
RsqCV = [RsqCV1 RsqCV2];
b = [beta1 beta2];

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [YHAT,MSE,Rsq] = getnnfits(Y,X,trainvec,testvec)
%%% think about returning *all* the YHAT, even though you only do the Rsq
%%% calculation on the held-out data


% fit a neural net
net = nnFit(X(trainvec,:)',Y(trainvec,:)');
YHAT = (sim(net,X(testvec,:)'))';

% get error stats
E = Y(testvec,:) - YHAT;
MSE = mean(E.^2);
MST = var(Y(testvec,:));
Rsq = 1 - MSE./MST;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function gainplots(gx,gy)

% setColors

for i = 1:2 % length(params.mods)
    figure()
    sinds = 1:1000;
    scatter(gx(sinds,i),gy(sinds,i))
    v = axis;
    a = max(v(1),v(3));
    b = min(v(2),v(4));
    hold on;
    plot(a:.1:b,a:.1:b,'k','LineWidth',2)
    hold off;
    xlabel('gain')
    ylabel('(normalized) input total spike count');
end



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [etaHAT,Rsq] = improveEta(p,q,trainvec,testvec)


% find the best linear fits
[etahatLR,~,~,RsqLR] = linearfits2(p.ttlSpks,q.ttlSpks,trainvec,testvec);

% find the best neural-net (nonlinear) fits
[etahatNN, ~, RsqNN] = getnnfits(p.ttlSpks,q.ttlSpks,trainvec,testvec);

% plot
figure; bar([RsqLR RsqNN]);
% ylabel('${\rm R}^2$','Interpreter','Latex','Fontsize',15)
ylabel('R^2','Interpreter','none');
ax = axis;  ax(3:4) = [0 1];    axis(ax);

% malloc
etaHAT = zeros(length(testvec),2);
Rsq = zeros(size(RsqLR));

% use the better fits
for j = 1:length(RsqLR)
    if RsqLR(j) > RsqNN(j),
        fprintf('using linear fit for eta(%d)\n',j);
        Rsq(j) = RsqLR(j);
        etaHAT(:,j) = etahatLR(:,j);
    else
        fprintf('using neural-network fit for eta(%d)\n',j);
        Rsq(j) = RsqNN(j);
        etaHAT(:,j) = etahatNN(:,j);
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function InfoINT = updatePrcns(p,ttlSpksHat,testvec)
%%% hard codes two modalities

InfoINT = sum(p.Info(testvec,:,:,1:2).*...
    permute(ttlSpksHat./p.ttlSpks(testvec,:),[1,3,4,2]),4);


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotPostCovs(varargin)
%%% make fontsize an argument??
%%% make number of figures to save an argument??

n = 4;

setColors;
clrs = [OPTcolor; VIScolor; PROPcolor; EFHcolor; EFHcolor];
stys = '----:';
xrange = 0;
yrange = 0;

% malloc
fighandles = zeros(n,1);

for i = 1:size(varargin{1}.Xpct,1)
    fh = figure; hold on;
    for j = 1:length(varargin)
        p = varargin{j};
        h = error_ellipse(inv(p(i).pcn),p(i).mu,'conf',.95);
        set(h,'LineWidth',1.5,'Color',clrs(j,:),'LineStyle',stys(j));
        xlabel('$\prop_1$ (rad)','Interpreter','Latex','Fontsize',25)
        ylabel('$\prop_2$ (rad)','Interpreter','Latex','Fontsize',25)
    end
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times-Roman','fontsize',25);
    
    if i <= n
        fighandles(i) = fh;
        axis('tight');
        ax = axis;
        xrange = max(xrange,ax(2)-ax(1));
        yrange = max(yrange,ax(4)-ax(3));
    end
    yorn = input('quit?','s');
    if strcmp(yorn,'y')
        break
    end
end

for i = 1:n
    set(0,'CurrentFigure',fighandles(i));
    ax = axis;
    ax(1) = ax(1) - (xrange - (ax(2)-ax(1)))/2;
    ax(2) = ax(2) + (xrange - (ax(2)-ax(1)))/2;
    ax(3) = ax(3) - (yrange - (ax(4)-ax(3)))/2;
    ax(4) = ax(4) + (yrange - (ax(4)-ax(3)))/2;
    axis(ax);
    filename = ['C:\#code\RBM\results\postcov\postex',num2str(i)];
    saveas(fighandles(i),filename,'epsc');
    fix_dottedline([filename,'.eps']);
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [pVIS,pPROP,pINT] = separatePosteriors(p)

pVIS.Xpct = p.Xpct(:,:,strcmp(p.srcs,'Hand-Position'));
pVIS.Info = p.Info(:,:,:,strcmp(p.srcs,'Hand-Position'));
pPROP.Xpct = p.Xpct(:,:,strcmp(p.srcs,'Joint-Angle'));
pPROP.Info = p.Info(:,:,:,strcmp(p.srcs,'Joint-Angle'));
pINT.Xpct = p.Xpct(:,:,strcmp(p.srcs,'opt'));
pINT.Info = p.Info(:,:,:,strcmp(p.srcs,'opt'));

end
%-------------------------------------------------------------------------%