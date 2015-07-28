function gainerrors(D0,x0,G0,wts,params)
% GAINERRORS    mean squared errors for the gain proxies
%   GAINERRORS computes and plots (&c.) the R^2s for the gain
%   proxies---etas for the input layer, some other thing for the hidden 
%   layer.
% clear; clc;
% load results/new/wtsBigSpace.mat
% load results/numhidswts/Std050.mat
% [D0 x0 G0] = DATAGENPP(1000,params);
%-------------------------------------------------------------------------%
% Revised: 10/01/12
%   -changed rhoZ to use (N-1) rather than N
% Revised: 01/30/12
%   -functionized etc.
%   -added function to compute covariance-distance errors
% Revised: 01/27/12
%   -rationalized (this function is still a mess)
% Created: 10/12/11 
%   by JGM
%-------------------------------------------------------------------------%

% init
params.smpls = 15;
[Ncases,Nvis,Nbatches] = size(D0);
Nexamples = Nbatches*Ncases;

% get gain proxies and plot
[g0, g1, g2, g3, d0, d1, d2, x0, rhoZ] = getgainproxies(D0,x0,G0,wts,params);
gainplots(g0,g1); gainplots(g1,g2);
% etahat = pseudodecode(v0,wts);

% find the best linear fits, g1->g0, g2->g1
[~, beta,MSElinear10,Rsq10] = linearfits(g0,g1);
[~, beta,MSElinear20,Rsq20] = linearfits(g0,g2);
[~, beta,MSElinear30,Rsq30] = linearfits(g0,g3);
[~, beta,MSElinear21,Rsq21] = linearfits(g1,g2);
[~, beta,MSElinear31,Rsq31] = linearfits(g1,g3);

% find the best neural-net (nonlinear) fits and plot
[g10hat,MSEnn10,Rsqnn10] = getnnfits(g0,g1);    % compare: eta(r) to gains
% [g21hat MSEnn21] = getnnfits(g1,g2);          %       eta(v) to eta(r)
[g20hat,MSEnn20,Rsqnn20] = getnnfits(g0,g2);    %       eta(v) to gains
[g30hat,MSEnn30,Rsqnn30] = getnnfits(g0,g3);    %       eta(vbar) to gains
bar([Rsqnn10 Rsqnn20 Rsqnn30]);


% get posterior covariances and their distances
trainvec = 1:Nexamples/2; testvec = trainvec + Nexamples/2;
% [covdist1 covdist2 covdist3]
[Rsq,Rsq2] = getpostcovs(g0(testvec,:),x0(testvec,:),rhoZ,params,...
    g10hat,d0(testvec,:),g20hat,d1(testvec,:),g30hat,d2(testvec,:));
%%%% may want to try both eta and eta+1??
5
%%% the unexplained thing that remains is that the input yields better
%%% gains, but essentially identically good covariances---and in fact much
%%% worse in the fourth component, which seems pretty dubious.

% scatter(g1hatnn(1,1:1500),g1(testvec(1:1500),1),[],VIScolor)


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [YHAT,b,MSE,Rsq] = linearfits(Y,X)

% init
N = size(Y,1);
trainvec = 1:N/2;
testvec = N/2 + trainvec;

% fit an affine ("linear") model
beta1 = regress(Y(trainvec,1),[X(trainvec,1) ones(N/2,1)]);
beta2 = regress(Y(trainvec,2),[X(trainvec,2) ones(N/2,1)]);
YHAT = [[X(testvec,1) ones(N/2,1)]*beta1,[X(testvec,2) ones(N/2,1)]*beta2];
b = [beta1 beta2];

% get error stats
E = Y(testvec,:) - YHAT;
MSE = mean(E.^2);
MST = var(Y(testvec,:));
Rsq = 1 - MSE./MST;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [YHAT,MSE,Rsq] = getnnfits(Y,X)

% init
N = size(Y,1);
trainvec = 1:N/2;
testvec = N/2 + trainvec;

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

for i = 1:2 % params.Nmods
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
function [MINSpostcov,xx] = postcovs(eta,dd,xx,params)

% init
Nexamples = size(eta,1);
DECODE = 0;

% malloc if necessary
if isempty(xx)
    xx = zeros(Nexamples,params.Ndims*params.Nmods);
    DECODE = 1;
else
    if ~isempty(dd)
        fprintf('warning: expected an empty arg; using d rather than x\n');
    end
end

% init
tuningCov = computetuningcovs(params);
Nmods = params.Nmods;
Ndims = params.Ndims;

% malloc
MINSpostcov = zeros(Ndims,Ndims,Nexamples);

% loop
for iExample = 1:Nexamples
    if DECODE
        shatL = decoder(dd(iExample,:),params);
        xx(iExample,:) = shatL(:)';
    end
    
    postprecision = 0;
    for iMod = 1:Nmods
        J = ntrlJacobian(xx(iExample,:),iMod,params);
        SILSpostcov = tuningCov{iMod}/eta(iMod);
        SINSpostcov = J*SILSpostcov*J';
        postprecision = postprecision + inv(SINSpostcov);
    end
    MINSpostcov(:,:,iExample) = inv(postprecision);
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [g0, g1, g2, g3, d0, d1, d2, x0, rhoZ] =...
    getgainproxies(D0,x0,G0,wts,params)

% init
Nmods = params.Nmods;
Ndims = params.Ndims;
N = params.N;
C = params.C;
gridsize0 = params.gridsize;
[Ncases,Nvis,Nbatches] = size(D0);
Nexamples = Nbatches*Ncases;
% numhids = params.numsUnits(2);


% malloc
% V0 = zeros(numcases,numhids,numbatches);
D1 = zeros(size(D0));
D2 = zeros(size(D0));
eta0 = zeros(Ncases,Nmods,Nbatches);
eta1 = zeros(Ncases,Nmods,Nbatches);
eta2 = zeros(Ncases,Nmods,Nbatches);

% get "updated inputs"
for i = 1:Nbatches
    means = feedforward(D0(:,:,i),wts{1}(1:end-1,:),wts{1}(end,:),...
        params.typeUnits{2},params);
    FXN = params.typeUnits{2};
    
    NNNN = params.smpls;
    states = sampler(means,FXN,params);      % samples!
    for j=1:NNNN-1
        states = states + sampler(means,FXN,params);
    end
    numspikes(i) = mean(sum(states,2));
    % numspikes = mean(mean(sum(V0,2),3));
    
    states = states/NNNN;
    % V0(:,:,i) = states;
    
    D1(:,:,i) = feedforward(states,wts{2}(1:end-1,:),wts{2}(end,:),...
        params.typeUnits{1},params);
    D2(:,:,i) = feedforward(means,wts{2}(1:end-1,:),wts{2}(end,:),...
        params.typeUnits{1},params);
    
    fprintf('.');
end


% total spike counts
for i = 1:Nmods
    eta0(:,i,:) = sum(D0(:,(1:N^Ndims) + (N^Ndims)*(i-1),:),2);
    eta1(:,i,:) = sum(D1(:,(1:N^Ndims) + (N^Ndims)*(i-1),:),2);
    eta2(:,i,:) = sum(D2(:,(1:N^Ndims) + (N^Ndims)*(i-1),:),2);
end


% The normalizing term that relates gains and etas; see PPCexpectedFI, e.g.
rho = ((N-1)/gridsize0)^Ndims;
rhoZ = rho*((2*pi)^(Ndims/2))*sqrt(det(C));

[g0, g1, g2, g3, d0, d1, d2, x0] =...
    longdata(G0,(eta0+1)/rhoZ,(eta1+1)/rhoZ,(eta2+1)/rhoZ, D0, D1, D2, x0);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function etahat = pseudodecode(v0,wts)
% use the *right* pseudoinverse (confusing b/c logit(v)-b = W'*u, w/' op)
%   => u = W*(W'*W)*(logit(v)-b)
%   => eta = M*W*(W'*W)*(logit(v)-b)

% init
N = size(wts{1},1);
Nexamples = size(v0,1);

% get right pseudoinverse
W = wts{1}(1:end-1,:);
b = wts{1}(end,:);
M = [ones(1,N^2) zeros(1,N^2); zeros(1,N^2) ones(1,N^2)];
R = M*W/(W'*W);

% invert
logit = @(x)(log(x./(1-x)));
etahat = (logit(v0)-repmat(b,Nexamples,1))*R';


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function varargout = getcovdistances(M,varargin)

% init
Nexamples = size(M,3);
COLOR = 'rgbcy';

% malloc
for i = 1:nargin-1
    covdistance{i} = zeros(Nexamples,1);
end

% compute covariance distances
for i = 1:nargin-1
    for j = 1:Nexamples
        covdistance{i}(j) = covdist(M(:,:,j),varargin{i}(:,:,j));
    end
end

% plot
figure; hold on;
cvr = squeeze(mean(M,3));
error_ellipse(cvr,[0;0],'style','k');
for i = 1:nargin-1 
    cvr = squeeze(mean(varargin{i},3));
    error_ellipse(cvr,[0;0],'style',COLOR(i));
end
axis equal; hold off;

% store the output
varargout = covdistance;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Rsq Rsq2] = getpostcovs(gTrue,xTrue,rhoZ,params,varargin)

% compute posterior covariances
%%% use the mode of the Poisson dstrb???
%%% MINSpostcovTrue = postcovs(floor(gTrue*rhoZ),[],xTrue,params);
[MINSpostcovTrue, ~] = postcovs(gTrue*rhoZ-1,[],xTrue,params);



% init
n = (nargin - 4)/2;
[dim1 dim2 N] = size(MINSpostcovTrue);

% malloc
Rsq = zeros(n,dim1*dim2);
Rsq2 = zeros(n,1);
% d = zeros(N,1);
% MINSpostcov = cell(1,n);
% covdst = cell(1,n);
KL = zeros(n,N);

% loop
PCparams0 = reshape(MINSpostcovTrue,dim1*dim2,N);
MST = var(PCparams0,[],2)';
clrs = 'rgb';
for i=1:n
    
    %%% use the mode of the Poisson dstrb???
    %%% MINSpostcov{i} = postcovs(floor(varargin{2*i-1}*rhoZ),varargin{2*i},[],params);
    % MINSpostcov{i} = postcovs(varargin{2*i-1}*rhoZ-1,varargin{2*i},[],params);
    % MINSpostcov{i} = postcovs(varargin{2*i-1}*rhoZ-1,[],xTrue,params);
    [MINSpostcov shat] = postcovs(varargin{2*i-1}*rhoZ-1,varargin{2*i},[],params);
    
    % compute entrywise R^2s while you're at it
    % PCparams1 = reshape(MINSpostcov{i},dim1*dim2,N);
    PCparams1 = reshape(MINSpostcov,dim1*dim2,N);
    err = (PCparams1 - PCparams0)';
    MSE = mean((err).^2);    
    Rsq(i,:) = 1 - MSE./MST;
    
    % compute KL divergence!
    for j = 1:N
        %%% hard-coded 3:4---fix
        KL(i,j) = gaussKLdvrg(xTrue(j,3:4)',xTrue(j,3:4)',... % shat(j,3:4)',...
            MINSpostcovTrue(:,:,j),MINSpostcov(:,:,j));
        
        if j < 0
            f2 = figure(j);
            set(0,'CurrentFigure',f2);
            hold on;
            h = error_ellipse(MINSpostcov(:,:,j),shat(j,3:4),'conf',0.95);
            set(h,'LineWidth',1,'Color',clrs(i));
            hold off;
        end
    end
end


%%% idea: measure avg. KL from the two unimodal's to the optimal posterior, 
%%% then do the same for the RBM-based posterior.  Also do this w/the
%%% unimodal variance but correct posterior mode, b/c the wrong posterior
%%% mode appears to contribute most of the distance.
%%% (2) why is the mean-based one so far away??


% % compute the covariance distances (see Mitteroecker 2009)
% [covdst{:}] = getcovdistances(MINSpostcovTrue,MINSpostcov{:});
% 
% % get R-squared-like measure
% for i = 1:size(MINSpostcovTrue,3)
%     d(i) = covdist(MINSpostcovTrue(:,:,1),MINSpostcovTrue(:,:,i));
% end
% MST = var(d);
% for i = 1:length(covdst)
%     MSE = mean(covdst{i}.^2);
%     Rsq2(i) = 1 - MSE/MST;
% end

% varargout = covdst;

end
%-------------------------------------------------------------------------%
