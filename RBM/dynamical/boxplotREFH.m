function boxplotREFH(datatype,Ndims,VARIANTS)
% boxplotREFH   Box plot of MSEs for rEFH and related EM models
%
% USAGE:
%   boxplotREFH('LTI-PPC',1,{});
%
%   boxplotREFH('LTI-PPC',1,{'withEC'});
%
% This should recreate Figure 4 in Makin/Dichter/Sabes PLoS C.B. 2015, as
% well as Figure 4 in REFH_JMLR [under review].
%
% This function uses the 'VaryingHids' .mats.

%%% TO DO 
% (1) EM1, EM2 -> \KFone,\KFtwo
% (2) it would be nice if you had 20 "experiments" of *both* MODELs
%%% 

%-------------------------------------------------------------------------%
% Created: 08/18/15
%   -by JGM
%-------------------------------------------------------------------------%

% params
Xprmt0 = 1;
Nxprmt = 20; %%% hard-coded (which is bad)
xvec = Xprmt0:Nxprmt;
setColors
variantstr = arrayfun(@(i)(['_',VARIANTS{i}]),1:length(VARIANTS),...
    'UniformOutput',false);
MODEL = [sprintf('%iD_%s',Ndims,datatype)];


% plot number 1
h = figure(5542); clf;
hold on;



[MSEsREFH,Nstates,Nsensory] = plotREFHruns(MODEL,variantstr,xvec,EFHcolor);
MSEsEMbad = plotEMruns(MODEL,variantstr,xvec,Nstates-1,EMcolor);
MSEsEMgood = plotEMruns(MODEL,variantstr,xvec,Nstates,EMBESTcolor);
hold off;
fileprefix = [MODEL,variantstr{:}];
matlab2tikzWrapper([fileprefix,'_boxplotsA_',date],h);

% plot number 2
altBoxplot(MSEsREFH,MSEsEMbad,MSEsEMgood,MODEL,variantstr,xvec,Nstates,Nsensory);
%%% You may still want to edit the tikz pics:
% (1) just *before* end{axis}, add:
%   \draw [thick,decoration={brace,mirror},decorate] 
%   (axis cs:1,0.00013) -- node[below=7pt] 
%   {rEFH \# hiddens} (axis cs:12,0.00013);
% (2) to the "begin{axis}," add:
%   x tick label style={rotate=90},
%   clip=false,
% and comment out:
%   xlabel={rEFH \# hiddens}
% (3) fix the colors (mycolor -> EFHcolor, etc. etc.)


% plot number 3
if strcmp(MODEL,'1D_LTI-PPC')&&(Nstates==2)
    % the box plot in REFH_JMLR, comparing with the TRBM and RTRBM
    [MSEOPT,MSESENSORY] = getOPTandSENSORYMSEs(MODEL,variantstr);
    altEFHBoxplot(MSEOPT,MSESENSORY,MSEsEMbad,MSEsEMgood,MSEsREFH,...
        MODEL,variantstr,xvec,Nstates,Nsensory)
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [MSEs,Nstates,Nsensory] = plotREFHruns(MODEL,variantstr,xvec,yrclr)

% load the appropriate file
fileprefix = [getdir('data'),'RBMish/EFHs/wts_'];
algorithm = 'rEFH';
filesuffix = [variantstr{:},'_VaryingHids.mat'];
load([fileprefix,algorithm,'_',MODEL,filesuffix]);
MSEs = squeeze(ydataTensor(end,xvec,:))';

% Ns
Nsensory = length(params.mods)*params.N^params.Ndims;
Nruns = size(ydataTensor,3);

% plot
clf; hold on;
fprintf(['\n%s results for %i runs, %i experiments\n'],MODEL,Nruns,length(xvec));
bh = boxplot(MSEs,'positions',xvec*Nsensory,...
    'colors',yrclr,'boxstyle','outline','symbol','+');
for i=1:size(bh,1),for j=1:size(bh,2),set(bh(i,j),'linewidth',1.5); end;end
Nstates = size(params.dynamics.A,2);

xlabel('\# hiddens')
ylabel('mean square error (rad$^2$)');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function MSEs = plotEMruns(MODEL,variantstr,xvec,Nstates,yrclr)

% params
whisker = 1.5;
qrtilefunc = @(y)(prctile(y,[75,50])-prctile(y,[50,25]));
Nxprmt = length(xvec);

% load the appropriate file
fileprefix = [getdir('data'),'RBMish/EMparams/LDSOrd',num2str(Nstates),'_'];
filesuffix = [variantstr{:},'_VaryingHids.mat'];
load([fileprefix,MODEL,filesuffix]);
Nruns = length(allparams);
Nsensory = params.N^params.Ndims;
fprintf(['\nEM',num2str(Nstates),' results for %i runs\n'],Nruns);

% the canonical test data
fileprefix = [getdir('data'),'RBMish/testdata/data_'];
filesuffix = [variantstr{:},'.mat'];
load([fileprefix,MODEL,filesuffix],'Rtest','Xtest','Qtest');

% convert to cumulants
[ShatTest,ttlSpksTest] = decodeDataPPC(Rtest,Xtest,Qtest,params);
InfoTest = GTPNposteriorInfo(ttlSpksTest,params);
cmlntsTest = cumulantNeutralize(ShatTest,InfoTest,params);

% collect MSEs
MSEs = zeros(Nruns,1);
S = latents2stims(Xtest,Qtest.latent2stim,params.mods,params.Ndims);
NSind = strcmp(params.mods,params.NS);
for iRun = 1:Nruns
    %%%%name = ['EM$^',num2str(Nstates),'$'];
    pEM = KFposteriorization(cmlntsTest,Qtest,allparams(iRun),params);
    err = pEM.Xpct(:,:,NSind) - S(:,:,NSind);
    MSEs(iRun,1) = mean(err.^2);
end

% get quartile/outlier stats for plotting
Q = prctile(MSEs,[25,50,75]);
upperW = Q(3)+whisker*(Q(3)-Q(1));
lowerW = Q(1)-whisker*(Q(3)-Q(1));
outliers = MSEs((MSEs>upperW)|(MSEs<lowerW));
shadedErrorBar(xvec*Nsensory,repmat(MSEs,[1,Nxprmt]),...
    {@median,qrtilefunc},{'-k','Color',yrclr});
%%% scatter(zeros(size(outliers)),outliers,'+','MarkerEdgeColor',EYEcolor);
if ~isempty(outliers)
    plot(xvec*Nsensory,outliers*ones(1,Nxprmt),'--','Color',yrclr);
end
ax = axis;
axis([-1,ax(2),min([outliers;ax(3)]),max([outliers; ax(4)])]);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function altBoxplot(MSEsREFH,MSEsEMbad,MSEsEMgood,MODEL,variantstr,xvec,...
    Nstates,Nsensory)

Xprmt0 = xvec(1);
setColors

h = figure(5543); clf;
yrlabels{1} = ['EM$^',num2str(Nstates-1),'$'];
yrlabels{2} = ['EM$^',num2str(Nstates),'$'];
for i=1:length(xvec),yrlabels{i+2}=num2str(Nsensory*xvec(i)); end
bh = boxplot(cat(2,MSEsEMbad,MSEsEMgood,MSEsREFH),...
    'positions',[Xprmt0-2,Xprmt0-1,xvec],...
    'colors',[EMcolor;EMBESTcolor;repmat(EFHcolor,[length(xvec),1])],...
    'labels',yrlabels,'boxstyle','outline','symbol','+');
for i=1:size(bh,1),for j=1:size(bh,2),set(bh(i,j),'linewidth',1.5); end;end
hold on;
plot([Xprmt0-2,Xprmt0-1,xvec],min(MSEsEMbad)*ones(length(xvec)+2,1),...
    'Color',EMcolor)
plot([Xprmt0-1,xvec],min(MSEsEMgood)*ones(length(xvec)+1,1),...
    'Color',EMBESTcolor)
hold off;
xlabel('\# hiddens')
ylabel('mean square error (rad$^2$)');


% write
fileprefix = [MODEL,variantstr{:}];
matlab2tikzWrapper([fileprefix,'_boxplotsB_',date],h);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [MSEOPT,MSESENSORY] = getOPTandSENSORYMSEs(MODEL,variantstr)

fileprefix = getdir('data');
filesuffix = [variantstr{:},'_VaryingHids.mat'];
load([fileprefix,'RBMish/EFHs/wts_','rEFH_',MODEL,filesuffix],'params');

% the canonical test data
fileprefix = [getdir('data'),'RBMish/testdata/data_'];
filesuffix = [variantstr{:},'.mat'];
load([fileprefix,MODEL,filesuffix],'Rtest','Xtest','Qtest');

% convert to cumulants
[ShatTest,ttlSpksTest] = decodeDataPPC(Rtest,Xtest,Qtest,params);
InfoTest = GTPNposteriorInfo(ttlSpksTest,params);
pSENSORY = cumulantNeutralize(ShatTest,InfoTest,params);
pSENSORY.name = 'sensory';

% Kalman filter
LDSparamsTrue = getLDSparams(params.dynamics);
pOPT = KFposteriorization(pSENSORY,Qtest,LDSparamsTrue,params);

% compute errors
S = latents2stims(Xtest,Qtest.latent2stim,params.mods,params.Ndims);
NSind = strcmp(params.mods,params.NS);
err = pOPT.Xpct(:,:,NSind) - S(:,:,NSind);
MSEOPT = mean(err.^2);
err = pSENSORY.Xpct(:,:,NSind) - S(:,:,NSind);
MSESENSORY = mean(err.^2);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function altEFHBoxplot(MSEOPT,MSESENSORY,MSEsEMbad,MSEsEMgood,MSEsREFH,...
    MODEL,variantstr,xvec,Nstates,Nsensory)

Nxprmts = length(xvec);

% load the TRBM and RTRBM results
fileprefix = [getdir('data'),'RBMish/EFHs/wts_'];

algorithm = 'TRBM';
filesuffix = [variantstr{:},'_VaryingHids.mat'];
load([fileprefix,algorithm,'_',MODEL,filesuffix]);
MSEsTRBM = squeeze(ydataTensor(end,xvec,:))';

algorithm = 'RTRBM';
filesuffix = [variantstr{:},'_VaryingHids.mat'];
load([fileprefix,algorithm,'_',MODEL,filesuffix]);
MSEsRTRBM = squeeze(ydataTensor(end,xvec,:))';


% put the MSEs into cells
MSEsTRBMcell = mat2cell(MSEsTRBM,size(MSEsTRBM,1),ones(size(MSEsTRBM,2),1));
MSEsRTRBMcell = mat2cell(MSEsRTRBM,size(MSEsRTRBM,1),ones(size(MSEsRTRBM,2),1));
MSEsREFHcell = mat2cell(MSEsREFH,size(MSEsREFH,1),ones(size(MSEsREFH,2),1));
allMSEs = [{MSEsEMbad},{MSEsEMgood},MSEsTRBMcell,MSEsRTRBMcell,MSEsREFHcell];

% put the colors into cells
yrcolors{1} = 'EMclr';
yrcolors{2} = 'EMBESTclr';
for iXprmt = 1:Nxprmts, yrcolors = [yrcolors,{'TRBMclr'}]; end
for iXprmt = 1:Nxprmts, yrcolors = [yrcolors,{'RTRBMclr'}]; end
for iXprmt = 1:Nxprmts, yrcolors = [yrcolors,{'rbmclr'}]; end

% write the x tick labels
xticklabels{1} = '\KFone';
xticklabels{2} = '\KFtwo';
for iNumhid = 1:Nxprmts
    xticklabels = [xticklabels,{int2str(xvec(iNumhid)*Nsensory)}];
end


% extra plot bits
yrunderbrace = ['\draw [thick,decoration={brace,mirror},decorate]',...
    '(axis cs:3,0.00001) -- node[below=7pt]',...
    '{rEFH/TRBM/RTRBM \# hiddens} (axis cs:22,0.00001);%'];
SENSORYplot = ['\addplot[color=pclr,line width=1.5,domain=1:22] {',...
    num2str(MSESENSORY),'};\addlegendentry{\KFnaught}%'];
% EMplot = ['\addplot[color=EMclr,line width=1.5,domain=1:22] {',...
%     num2str(min(MSEsEMbad)),'};%'];
EMplot = ['\addlegendimage{no markers, EMclr, line width=1.5}',...
    '\addlegendentry{\KFone}%'];
% EMBESTplot = ['\addplot[color=EMBESTclr,line width=1.5,domain=2:22] {',...
%    num2str(min(MSEsEMgood)),'};%'];
EMBESTplot = ['\addlegendimage{no markers, EMBESTclr, line width=1.5}',...
    '\addlegendentry{\KFtwo}%'];
OPTplot = ['\addplot[color=optclr,line width=1.5,domain=1:22] {',...
    num2str(MSEOPT),'};\addlegendentry{\KFopt}%'];
TRBMplot =['\addlegendimage{no markers, TRBMclr, line width=1.5}',...
     '\addlegendentry{TRBM}%'];
RTRBMplot = ['\addlegendimage{no markers, RTRBMclr, line width=1.5}',...
    '\addlegendentry{RTRBM}%'];
rEFHplot = ['\addlegendimage{no markers, rbmclr, line width=1.5}',...
    '\addlegendentry{rEFH}%]'];


% extra axis options
rotatexlabels = 'x tick label style={rotate=90}';
setwidth = 'width = 6.0in';
setheight = 'height = 4.0in';
allMSES = [MSEOPT;MSESENSORY;MSEsEMbad;MSEsEMgood;...
    MSEsREFH(:);MSEsTRBM(:);MSEsRTRBM(:)];
MSEmin = min(allMSES); MSEmax = max(allMSES);
MSErange = MSEmax-MSEmin;
MSEmin = MSEmin - 0.02*MSErange;
MSEmax = MSEmax + 0.02*MSErange;
setdims = ['xmin=0, xmax=23, ymin=',num2str(MSEmin),', ymax=',num2str(MSEmax)];


% write out a tikz/pgf plot
tikzBoxPlot([1,2,2+(1:Nxprmts),2+(1:Nxprmts),2+(1:Nxprmts)],allMSEs,...
    yrcolors,xticklabels,'mean square error (rad$^2$)',...
    {rotatexlabels,setwidth,setheight,setdims},...
    {SENSORYplot,EMplot,EMBESTplot,OPTplot,TRBMplot,RTRBMplot,rEFHplot,...
    yrunderbrace},...
    [MODEL,'_ErrorStatsVsNumHiddens_',date]);




end
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function altEFHBoxplotOld(MSEOPT,MSESENSORY,MSEsEMbad,MSEsEMgood,MSEsREFH,MODEL,...
    xvec,Nstates,Nsensory)

Xprmt0 = xvec(1);

% load the TRBM and RTRBM results
fileprefix = [getdir('data'),'RBMish/EFHs/wts_'];
filesuffix = '_VaryingHids.mat';
load([fileprefix,'1DTRBM',filesuffix]);
MSEsTRBM = squeeze(ydataTensor(end,xvec,:))';
load([fileprefix,'1DRTRBM',filesuffix]);
MSEsRTRBM = squeeze(ydataTensor(end,xvec,:))';

% init colors
setColors
% for rEFH.pdf, colors are changed....
TRBMcolor = [190,186,218]/255;                % purple
RTRBMcolor = [251,128,114]/255;               % red-ish
EMcolor = [179,222,105]/255;                  % yellow
EMBESTcolor = [128,177,211]/255;              % blue


% begin figure
h = figure(8001); clf; hold on;
yrlabels{1} = ['EM$^',num2str(Nstates-1),'$'];
yrlabels{2} = ['EM$^',num2str(Nstates),'$'];
for i=1:length(xvec),yrlabels{i+2}=num2str(Nsensory*xvec(i)); end


% TRBM
bh = boxplot(cat(2,MSEsTRBM),...
    'positions',xvec,...
    'colors',repmat(TRBMcolor,[length(xvec),1]),...
    'labels',yrlabels(3:end),'boxstyle','outline','symbol','+');
for i=1:size(bh,1),for j=1:size(bh,2),set(bh(i,j),'linewidth',1.5); end;end

% RTRBM
bh = boxplot(cat(2,MSEsRTRBM),...
    'positions',xvec,...
    'colors',repmat(RTRBMcolor,[length(xvec),1]),...
    'labels',yrlabels(3:end),'boxstyle','outline','symbol','+');
for i=1:size(bh,1),for j=1:size(bh,2),set(bh(i,j),'linewidth',1.5); end;end

% EMbad, EMgood, and rEFH
bh = boxplot(cat(2,MSEsEMbad,MSEsEMgood,MSEsREFH),...
    'positions',[Xprmt0-2,Xprmt0-1,xvec],...
    'colors',[EMcolor;EMBESTcolor;repmat(EFHcolor,[length(xvec),1])],...
    'labels',yrlabels,'boxstyle','outline','symbol','+');
for i=1:size(bh,1),for j=1:size(bh,2),set(bh(i,j),'linewidth',1.5); end;end

% *lines* for becnhmarks (SENSORY, EMbad, EMgood, OPT)
plot([Xprmt0-2,xvec(end)],MSEOPT*ones(2,1),'Color',OPTcolor);
plot([Xprmt0-2,xvec(end)],MSESENSORY*ones(2,1),'Color',PROPcolor);
plot([Xprmt0-2,xvec(end)],min(MSEsEMbad)*ones(2,1),'Color',EMcolor);
plot([Xprmt0-1,xvec(end)],min(MSEsEMgood)*ones(2,1),'Color',EMBESTcolor);

% fix the axes
allMSES = [MSEOPT;MSESENSORY;MSEsEMbad;MSEsEMgood;...
    MSEsREFH(:);MSEsTRBM(:);MSEsRTRBM(:)];
MSEmin = min(allMSES);
MSEmax = max(allMSES);
MSErange = MSEmax-MSEmin;
MSEmax = MSEmax + 0.02*MSErange;
MSEmin = MSEmin - 0.02*MSErange;
ax = axis;
axis([ax(1:2),MSEmin,MSEmax])

% label and save
xlabel('\# hiddens')
ylabel('mean square error (rad$^2$)');
hold off;
matlab2tikzWrapper([MODEL,'boxplotsB',date],h);


end
%-------------------------------------------------------------------------%




















