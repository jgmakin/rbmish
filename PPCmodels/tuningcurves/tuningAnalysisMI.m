function tuningAnalysisMI(wts,params)
% tuningAnalysisMI    Analyzes the standard network's hidden units
%   tuningAnalysisMI(wts,params) analyzes the tuning of the hidden neurons 
%   of the standard (multisensory-integration) network, including showing
%   tuning curves for the 1D manifold over which LMMM (McGuire2011) trained
%   her monkey.
%
%   USAGE:
%   load results/numhidswts/Std050.mat
%   tuningAnalysisMI(wts,params);


%-------------------------------------------------------------------------%
% Revised: 01/01/17
%   -changed to work with *functions* stored in params rather than
%   generateData.m.
% Revised: 04/01/13
%   -cleaned
% Revised: 01/24/13
%   -rationalized, functionized, renamed from LMMMtargs
% Created: 01/xx/13
%   by JGM
%-------------------------------------------------------------------------%

%%%%%%%%%%%
% NB: this runs, but you haven't checked all the figures
%   -jgm, 1/1/17
%%%%%%%%%%%


% init params
gmin = params.gmin;
gmax = params.gmax;
params.gmin = mean([gmin; gmax]);
params.gmax = mean([gmin; gmax]);

% generate data on an arc
[S,phi] = generateArcData(params);

% get the workspace and plot it
wrkspc = getWorkspace(100,params);
plotArcAndWorkspace(S,wrkspc,length(params.mods));

% get wts (and params)
wts0 = initializeWts(params);


% generate noiseless V given s and g
V12 = noiselessVgivenSandG(S,[12 12],wts,params);    % tuned
V15 = noiselessVgivenSandG(S,[15 15],wts,params);
V18 = noiselessVgivenSandG(S,[18 18],wts,params);

V012 = noiselessVgivenSandG(S,[12 12],wts0,params);  % untuned
V015 = noiselessVgivenSandG(S,[15 15],wts0,params);
V018 = noiselessVgivenSandG(S,[18 18],wts0,params);


% find the corresponding tuning curve for the whole space
S = getStimuliTiled(100^2,class(V12),params);
Vall = noiselessVgivenSandG(S,[15 15],wts,params);
V0all = noiselessVgivenSandG(S,[15 15],wts0,params);



% count the number of neurons on under these different conditions; maybe
% average across stimuli.
[avgNumberAboveChance12,avgNumberON12] = percentOn(S,12,15,wts,params);
[avgNumberAboveChance15,avgNumberON15] = percentOn(S,15,15,wts,params);
[avgNumberAboveChance18,avgNumberON18] = percentOn(S,18,15,wts,params);



% plot
tunedneurons = plotHiddenTuningCurves({V12,V15,V18},{V012,V015,V018},...
    Vall,V0all,phi,S(:,:,strcmp(params.mods,'Joint-Angle')),params);
plotIncomingWts(tunedneurons,wrkspc,wts,params);





end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [S,phi] = generateArcData(params)

phi = linspace(120,60,params.Ncases*15);
xx = 26*[cosd(phi)', sind(phi)'] + repmat([-5,0],[length(phi),1]);
th = IK2link(xx,params.roboparams,1);
S = cat(3,xx,th);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotArcAndWorkspace(ArcStims,wrkspc,Nmods)

for iMod = 1:Nmods
    figure(iMod); hold on;
    plot(wrkspc{iMod}(:,1),wrkspc{iMod}(:,2),'k','LineWidth',1.5);
    scatter(ArcStims(:,1,iMod),ArcStims(:,2,iMod),'k.');
    hold off;
end

end
%-------------------------------------------------------------------------%

    
%-------------------------------------------------------------------------%
function wts = initializeWts(params)

Nvis = sum(params.numsUnits{1});
Nhid = sum(params.numsUnits{2});

vishid      = 0.01*randn(Nvis,Nhid);
visbiases   = 0*ones(Nvis,1);
hidbiases   = 0*ones(1,Nhid);

wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases'];


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function V = noiselessVgivenSandG(S,g,wts,params)

% init
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};
[Nexamples,Ndims,Nmods] = size(S);

Q.G = repmat(g,[Nexamples,1]);
Q.biases = zeros([Ndims,Nmods],'like',S);
Q.DECOUPLED = zeros(Nexamples,1,'like',S);
R = params.getData(S,Q);
%%% R = generateData(size(S,1)/params.Ncases,params,'stimuli',S);
V = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,hidNums,params);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function tunedneurons = plotHiddenTuningCurves(V,V0,Vall,V0all,phi,th,params)

% init
n = 4;
Nhid = sum(params.numsUnits{2});

% plot the first n^2 hidden units
% figure(3); figure(4);
k = 0;
% vv = shortdata(M,V2);
tunedneurons = [];


figure(3); clf; p = panel(); p.pack(n,n);
figure(7); clf; q = panel(); q.pack(n,n);
figure(4); clf; r = panel(); r.pack(n,n);
figure(8); clf; s = panel(); s.pack(n,n);
for i = 1:Nhid
    if (max(V{1}(:,i))>0.05)&&(min(V{1}(:,i)<0.95))&&(k<n^2)
        % 1 b/c we're only checking g = 12, the worst case
        
        % this is a tuned neuron
        tunedneurons = [tunedneurons; i];
        k = k+1;
        
        % plot tuning curves along the arc of phi
        plotArcTuningCurves(p,n,k,i,3,phi,V,[]);
        plotArcTuningCurves(q,n,k,i,7,phi,V0,[]);
        
        % plot tuning curves for the whole space
        plotEntireTuningCurves(r,n,k,i,4,th,Vall,params,[]);
        plotEntireTuningCurves(s,n,k,i,8,th,V0all,params,[]);
        
    end
end
p.de.margin = 2;
q.de.margin = 2;
r.de.margin = 2;
s.de.margin = 2;


% tunedneurons = [27 32 43 47 59 62 92 117 129 130 136 143,...
%     145 149 166 169 172 179 188 189 198 227 234 245 258]


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotIncomingWts(tunedneurons,wrkspc,wts,params)

% init
n = 4;
smin = params.smin;
smax = params.smax;
N = params.N;

% prepare the boundaries
%%% malloc
for j = 1:2
    xi(:,j) = grid2world([1 1],[smin(:,j),smax(:,j)],params);
    xf(:,j) = grid2world([N N],[smin(:,j),smax(:,j)],params);
end

figure(5); clf; p{1} = panel(); p{1}.pack(n,n);
figure(6); clf; p{2} = panel(); p{2}.pack(n,n);
for i = 1:length(tunedneurons)
    
    % plot incoming wts
    W = displayshape(wts{2}(tunedneurons(i),:),params);
    
    % for both modalities
    for j = 1:2
        figure(4+j); 
        % subplot(n,n,i);
        p{j}(floor((i-1)/n)+1,mod(i-1,n)+1).select();
        imagesc([xi(1,j);xf(1,j)],[xi(2,j);xf(2,j)],W{j}');
        hold on;
        plot(wrkspc{j}(:,1),wrkspc{j}(:,2),'k','LineWidth',1.5);
        hold off;
        axis xy
        colormap(gray)
        if (floor((i-1)/n)+1)~=n, set(gca, 'xtick', []); end
        if (mod(i-1,n)+1)~=1, set(gca, 'ytick', []); end
        axis tight
    end
    
    
end
p{1}.de.margin = 2;
p{2}.de.margin = 2;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotArcTuningCurves(p,n,k,i,fignum,phi,V,titlestring)

clrs = 'bgrk';
figure(fignum);
% subplot(n,n,k); 
p(floor((k-1)/n)+1,mod(k-1,n)+1).select();
hold on;
for j = 1:length(V)
    try
    plot(90-phi,V{j}(:,i),clrs(j),'LineWidth',3);
    catch 
        keyboard
    end
    % scatter(90-phi,V{j}(:,i),[clrs(j),'.']);
end
axis([90-phi(1),90-phi(end),0,1])
if (floor((k-1)/n)+1)~=n, set(gca, 'xtick', []); end
if (mod(k-1,n)+1)~=1, set(gca, 'ytick', []); end
title(titlestring)
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plotEntireTuningCurves(p,n,k,i,fignum,th,Vall,params,titlestring)

% init
thmin = params.roboparams.thmin;
thmax = params.roboparams.thmax;
M = sqrt(size(Vall,1));

% plot tuning curves for the whole space
figure(fignum); 
% subplot(n,n,k);
p(floor((k-1)/n)+1,mod(k-1,n)+1).select();
hold on;
imagesc([thmin(1),thmax(1)],[thmin(2),thmax(2)],...
    reshape(Vall(:,i),M,M),[0,1]);
    % squeeze(vv(:,i,:)),[0 1]);  %%% this is apparently wrong
axis xy;
plot(th(:,1),th(:,2),'k.')
% scatter(th(1,:),th(2,:),8,'k.')
axis tight;
if (floor((k-1)/n)+1)~=n, set(gca, 'xtick', []); end
if (mod(k-1,n)+1)~=1, set(gca, 'ytick', []); end
title(titlestring)
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [avgNumberAboveChance,avgNumberON] =...
    percentOn(S,g,nSamples,wts,params)

V = noiselessVgivenSandG(S,[g g],wts,params);

avgNumberAboveChance = mean(mean(V > 0.5,2));

spikes = (repmat(V,[1 1 nSamples]) > rand([size(V),nSamples]));
avgNumberON = mean(mean(sum(spikes,3) > 0,2));



end
%-------------------------------------------------------------------------%


