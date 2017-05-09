function tuningAnalysisCX(wts,params)
% tuningAnalysisCX    Analyzes "addition" hidden units
%   tuningAnalysisCX(wts,params) analyzes the tuning of the hidden neurons
%   of the "addition" network after the fashion of Bremner & Andersen 2012
%   (originally in Pesaran 2006, 2010)....
%
% NB!!! That M=4 should be used to compare with those papers, but that M=40
% should be used to plot the sample tuning curves, in order to get decent
% resolution.

%-------------------------------------------------------------------------%
% Revised: 01/25/13
%   -added fxn drawTGspace to do just that
% Renamed: 01/24/13
%   -added "CX" at the end of the file name to distinguish from "MI."
%   -added fxn plotSampleTuningCurves to do just that
% Revised: 08/22/12
%   -added histogramming function to match Bremner 2012, fig. 4 (?)
% Revised: 08/06/12
%   -functionalized, rationalized
% Created: 08/03/12
%   by JGM
%-------------------------------------------------------------------------%

% init
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
gmin = params.gmin;
gmax = params.gmax;
params.gmin = mean([gmin; gmax]);
params.gmax = mean([gmin; gmax]);
%%% Nmods = length(params.mods);
%%% g = params.g;
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};
M = 4; %40
spacefignum = 101;

% scatter random samples on G,T space
scatterParallelogram(spacefignum,params);
plotParallelogram(spacefignum+1,200,params);


% init
eyespacefractions = [1/3, 1/2, 2/3];
clrs = 'rgb';
for i = 1:length(eyespacefractions)
    
    tic
    % generate S, R, and V
    [S,Q] = getStimuliTiled(M^2,dataclass,params,'eyefraction',eyespacefractions(i));
    params.typeUnits{1} = {'Dirac'};
    R = params.getData(S,Q);
    V = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,hidNums,params);
    
    % plot
    % tunedneurons = plotSampleTuningCurves(X,V,M);
    drawTGspace(eyespacefractions(i),spacefignum,clrs(i),params); 
    drawTGspace(eyespacefractions(i),spacefignum+1,clrs(i),params); 

    % perform gradient and separability analyses
    [res,gradIsSignif,separIsSignif] = gradAndSeparAnalysis(shortdata(M,3,V));

    % plot the gradients, labeled as they are separable and significant.
    plotResultant(res,gradIsSignif,separIsSignif)
    histResultant(res,gradIsSignif,separIsSignif)
    toc
    
end

6




% plotTCs(X(1:M:end,1),V(1:M:end,:),params);
% plotTCs(X(1:M,1),V(1:M,:),params);
%%% a hack---to get plotTCs to do plot these 1D vars in 2D
% params.m = 2;
% %%%
% plotTCs([X(:,1),X(:,3)],V,params)
% params.m = 1;



% (3) fix separability analysis---two big sv's should => not separable
%   [UPDATE: A threshold of sv1 > 4*sv2 reduces the proportion of
%   "separable" cells from 65% (78/120) to 20% (24/120), and in particular
%   makes all of the T-G cells inseparable.  Those latter are likely, then,
%   to be retinotopic; and the others are pure gain or pure target (i.e.,
%   body centered).  Many of the pure gain and pure target cells are also
%   "inseparable" under these criteria, however, which is surprising....
%   These turn out to be





% v = shortdata(M,V);
% ang = mod(angle(res(:,1) + 1i*res(:,2)),pi);
%
% x = ((pi - ang) < pi/8) + (ang < pi/8);
% y = gradIsSignif.*(~separIsSignif);
% z = logical(y.*x);
% w = find(z);
%
% for i = 1:length(w)
%     figure(11);
%     vv = squeeze(v(:,w(i),:));
%     imagesc(vv),
%     ss = svd(vv-mean(vv(:)));
%     fprintf('singular values: %f, %f; ratio: %f\n',ss(1),ss(2),ss(2)/ss(1));
%
%     plotResultant(res(w(i),:),1,1)
%     pause();
%     close
% end




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function scatterParallelogram(fignum,params)

% generate random data
Nexamples = 4000;
[~, S] = generateData(Nexamples,params);

% redefine to match Bremner2012
Tbody = squeeze(S(:,:,1)-S(:,:,3));
Gaze = squeeze(-S(:,:,3));

% plot
figure(fignum)
scatter(Gaze,Tbody,20,'k.')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotParallelogram(fignum,M,params)

gazemin = -params.roboparams.eyemax;
gazemax = -params.roboparams.eyemin;
posmin = params.roboparams.posmin;
posmax = params.roboparams.posmax;

gazepts = linspace(gazemin,gazemax,M);
pospts = linspace(posmin,posmax,M);
tbody(1,:) = [gazepts, gazemax*ones(1,M), fliplr(gazepts), gazemin*ones(1,M)];
tbody(2,:) = [gazepts + posmin, gazemax + pospts,...
    posmax + fliplr(gazepts), gazemin + fliplr(pospts)];

% plot
figure(fignum);
plot(tbody(1,:),tbody(2,:),'k');
axis equal;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [res,gradIsSignif,separIsSignif] = gradAndSeparAnalysis(vv)

% init
Nhid = size(vv,2);

% malloc
gradIsSignif = zeros(Nhid,1);
SVIsSignif = zeros(Nhid,1);
separIsSignif = zeros(Nhid,1);
res = zeros(Nhid,2);


% open for parallel processing
[pool,HADBEENCLOSED] = parallelInit;


% loop
tic
for k = 1:Nhid
    
    % the matrix for this hidden unit
    v = squeeze(vv(:,k,:));
    
    % gradient analysis
    res(k,:) = gradAnalysis(v);
    
    % separability analysis
    % [Vstar sv] = maxSepar(v);
    sv = maxSV(v);
    
    % randomization test
    [gradIsSignif(k), SVIsSignif(k)] = RandomizationTest(v,norm(res(k,:)),sv);
    
    % plot resultant gradient
    % plotResultant(res(k,:),gradIsSignif(k));
    
    % plot the singular values
    separIsSignif(k) = plotSVs(gradIsSignif(k)*SVIsSignif(k),v);
    
end
toc
if HADBEENCLOSED, delete(pool); end

% "separability" is significant if the first sv is large *and* gradIsSignif
% separIsSignif = gradIsSignif.*SVIsSignif;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [gradIsSignif,SVIsSignif] = RandomizationTest(v,magv,sv)
%%%%
% This could probably be vectorized
%%%%


% init
P = 10000;
M = size(v,1);
alpha = 0.95;

% define significant function
sigfxn = @(fv,fw)(sum(fv > fw)/P > alpha);

% for lots of trials...
parfor j = 1:P
    
    % randomize v
    [y,ind] = sort(rand(1,numel(v)));
    w = reshape(v(ind),M,M);
    
    % separability analysis
    % [Vstar sw(j)] = maxSepar(w);
    sw(j) = maxSV(w);
    
    % gradient analysis
    res = gradAnalysis(w);
    magw(j) = norm(res);
end

% is the gradient significantly bigger than a random one?
gradIsSignif = sigfxn(magv,magw);

% is the true first singular value greater than 95% of the fake ones?
SVIsSignif = sigfxn(sv,sw);



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function res = gradAnalysis(v)


[vx,vy] = gradient(v);

mag = sqrt(vx(:).^2 + vy(:).^2);
ang = atan2(vy(:),vx(:));               % angle \in [-pi,pi]
ang = mod(ang,pi);                      %       \in [0,pi]
ang = 2*ang;                            %       \in [0,2*pi]
x = mag.*cos(ang);                      % convert back to Cartesian
y = mag.*sin(ang);                      %
res = [mean(x) mean(y)];                % compute the average gradient


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function s1 = maxSV(v)

% mean subtract, then find svd
v = v - mean(v(:));
s = svd(v);
s1 = s(1);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotResultant(res,gradIsSignif,separIsSignif)

% init
phi = 0:0.05:2*pi;
figure();


hold on;
plot(cos(phi)/5,sin(phi)/5,'g'); hold on;
%%% divided arbitrarily by 5 to make visible

% loop through all hidden units
for k = 1:length(gradIsSignif)
    
    % set the color
    if gradIsSignif(k),
        if separIsSignif(k), clr = 'k'; else clr = 'g'; end
    else
        clr = 'r';
    end
    % if thing(k), clr = 'k'; else clr = 'r'; end
    
    % plot
    plot([0 res(k,1)], [0 res(k,2)],clr);
    scatter(res(k,1),res(k,2),[clr,'x']);
    axis equal;
    % pause(0.2)
end
fprintf('black have significant gradients and are separable\n');
fprintf('green have significant gradients and are inseparable\n');
fprintf('red have insignificant gradients\n');
hold off;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function histResultant(res,gradIsSignif,separIsSignif)

figure();
edges = linspace(-pi,pi,16);
%%% 16 increments to match Bremner 2012


inds = gradIsSignif&~separIsSignif;
ResInseparable = atan2(res(inds,2),res(inds,1));
NI = hist(ResInseparable,edges);

inds = gradIsSignif&separIsSignif;
ResSeparable = atan2(res(inds,2),res(inds,1));
NS = histc(ResSeparable,edges);

% plot
bar(edges,[NI(:)';NS(:)']','stacked')
set(gca,'XTick',-pi:pi/2:pi)
set(gca,'XTickLabel',{'G','T-G','T','T+G','G'})
legend('inseparable','separable')


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function separIsSignif = plotSVs(separIsSignif,v)

% init
thresh = 1/4; % 1
if separIsSignif, clr = 'k'; else clr = 'r'; end;

% get singular values
v = v - mean(v(:));
s = svd(v);

% plot
if 0
    plot(s,clr);
    hold on;
    scatter([1 2],s(1:2),[clr,'x']);
    hold off;
    
    % hang on...
    pause;
end

if (s(2)/s(1) > thresh), separIsSignif = 0; end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function tunedneurons = plotHiddenTuningCurves(p,x,V,tunedneurons,clr)


for j = 1:4
    for i = 1:3
        p(i,j).select();
        hold on;
        plot(x,V(:,tunedneurons(i+(j-1)*3)),'Color',clr);
        hold off;
        ax = axis;
        axis([ax(1) ax(2) 0 1]);
        
        if i~=3, set(gca, 'xtick', []); end
        if j~=1, set(gca, 'ytick', []); end
        
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function alltunedneurons = plotSampleTuningCurves(X,V,M)

% these turn out to be useful examples of each kind of shift (see paper)
goodvector = [62 63 30 43 27 52 42 64 38 19 54 8];

% Tret = Tbody - Gaze
Tbody = X(:,1)-X(:,3);
Gaze = -X(:,3);
GazeDiff = Gaze(end-M+1) - Gaze(1);

% get tuned neurons (in the first set of gains)
alltunedneurons = [];
for i=1:size(V,2)
    if (max(V(1:M,i)) - min(V(1:M,i))) > 0.1
        alltunedneurons = [alltunedneurons; i];
    end
end
tunedneurons = alltunedneurons(goodvector);
figure; clf; p = panel(); p.pack(3,4);
plotHiddenTuningCurves(p,Tbody(1:M),V(1:M,:),tunedneurons,'r');
plotHiddenTuningCurves(p,Tbody((end-M+1):end),V((end-M+1):end,:),tunedneurons,'g');
plotHiddenTuningCurves(p,Tbody(1:M)+GazeDiff,V(1:M,:),tunedneurons,'b');
p.de.margin = 2;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function drawTGspace(eyespacefraction,fignum,clr,params)
%%% you probably need to double check this for different values of posmin,
%%% posmax, eyemin, eyemax.  E.g., you probably need tmin.....

% get the corners of the parallelogram
emin = params.roboparams.eyemin;
emax = params.roboparams.eyemax;
% tmin = emin + params.posmin;
tmax = emax + params.roboparams.posmax;

% find the corners of this box-inscribed-in-the-parallelogram
m = (emax - emin)/tmax;
bmin(1) = emin*eyespacefraction;
bmax(1) = emax*eyespacefraction;
bmin(2) = m*(bmax(1) - emax);
bmax(2) = m*(bmin(1) - emin);

% get the box itself
MM = 100;
box(1,:) = [linspace(bmin(1),bmax(1),MM),bmax(1)*ones(1,MM),...
    linspace(bmax(1),bmin(1),MM),bmin(1)*ones(1,MM)];
box(2,:) = [bmin(2)*ones(1,MM),linspace(bmin(2),bmax(2),MM),...
    bmax(2)*ones(1,MM),linspace(bmax(2),bmin(2),MM)];

% plot the box
figure(fignum); hold on;
plot(box(1,:),box(2,:),clr,'LineWidth',1.5);
hold off;
axis equal

end
%-------------------------------------------------------------------------%







