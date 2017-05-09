function SLMD
% (Re)produce figures from Makin2017:
% "Superior Limb-Movement Decoding from Cortex with a New, Unsupervised-
%   Learning Algorithm"
%
% This file should generate all the (non-schematic) figures for the paper.
%
% To reacquire the filters, including retraining all the rEFHs, run 
% BMImaster.  Alternatively, you can set TRAINNEWREFHS and/or TRAINNEWLDS 
% to zero and used saved filters, but re-acquire the R^2s
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Revised: 04/15/17
% Revised: 03/20/17
%   -functionized, rationalized
% Created: 02/04/17
%   -by JGM
%-------------------------------------------------------------------------%


%%%%
% TO DO: 
% (2) add \KFobs and \KFem via \providecommand at the top of the tikz file,
% and then change this file accordingly....
%%%%

% make sure that all figs get saved in the same, predictable dimensions
set(0,'DefaultFigureWindowStyle','normal')

% "metadata"
monkeys = {'Indy','Chewie','Loco'};

% this sets the names but also the "canonical" ordering!
[decoders(1:7).name] = deal('static','kfobs','kfemstatic',...
    'kfemdynamic','ukf','refhstatic','refhdynamic');

% load the saved results
[decoders,kinemat,binwidths,trainingtimes,NdataTest,Nneurons,mkinds] =...
    assembleData(decoders,monkeys,'Rsqs_allBins');
fprintf('"canonical" training time is %d seconds\n',trainingtimes);




keyboard


% rEFH vs. nearest rival
[B,p] = nearestRivalPlots(decoders,kinemat,binwidths,mkinds,monkeys,...
    'refhdynamic',{'ukf'});

% how much improvement for pos, vel, acc?
summaryStats(decoders,kinemat,binwidths,NdataTest,'refhdynamic',{'kfobs','ukf'});

% bar summary plots (\FigAllDecodersBarPlots)
allDecodersBarPlots(decoders,kinemat,binwidths,NdataTest);

% effect of binwidths
effectOfBinwidthPlots(decoders,kinemat,binwidths,NdataTest,0);

% effect of number of training samples
effectOfNumberOfTrainingSamples(decoders,monkeys,0);

% effect of dropping out neurons
effectOfNumberOfNeurons(decoders,kinemat,binwidths,Nneurons,mkinds,monkeys,0);

% effect of rEFH components
effectOfREFHcomponents(monkeys);







end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [decoders,kinemat,binwidths,trainingtimes,NdataTest_all,...
    Nneurons_all,mkinds,varargout] = assembleData(decoders,monkeys,fileprefix)
% R^2 and SNR: rEFH vs. nearest rival (\FigNearestRivalScatterPlots)


% init
NdataTest_all   = [];
Nsessions       = 0;
Nneurons_all    = [];
allRsqs         = [];

% for each monkey...
for iMk = 1:length(monkeys)
    
    % load its data
    load(sprintf('%s/RBMish/BMI/%s_%s.mat',...
        getdir('data'),fileprefix,monkeys{iMk}));
    if iMk > 1
        if ~all(strcmp({kinemat.name},{kinemat_prev.name}))
            error('mismatching kinematic variables in monkey %i',iMk)
        end
        if ~all(binwidths == binwidths_prev)
            error('mismatching binwidths in monkey %i',iMk)
        end
        if ~all(trainingtimes == trainingtimes_prev)
            error('mismatching training times in monkey %i',iMk)
        end
        if exist('fracneurons','var')
            if ~all(fracneurons == fracneurons_prev)
                error('mismatching fractions of neurons in monkey %i',iMk)
            end
        end
    end
    
    % for comparison the next time through
    kinemat_prev = kinemat;
    binwidths_prev = binwidths;
    trainingtimes_prev = trainingtimes;
    if exist('fracneurons','var'), fracneurons_prev = fracneurons; end
    
   
    % accumulate Rsqs
    [~,iDecoders] = ismember({decoders(:).name},decodernames);
    allRsqs = cat(1,allRsqs,Rsqs(:,:,:,iDecoders));
    
    % accumulate other useful numbers
    Msessions       = size(Rsqs,1);
    mkinds{iMk}     = Nsessions + (1:Msessions);
    Nsessions       = Nsessions + Msessions;
    NdataTest_all   = [NdataTest_all; NdataTest];
    if exist('Nneurons','var'), Nneurons_all=[Nneurons_all;Nneurons]; end
    
end

% now break apart the big tensor to assign CoDs to each decoder
for iDecoder = 1:length(decoders)
    decoders(iDecoder).CoD = allRsqs(:,:,:,iDecoder);
end

% also assign the other useful bits for these decoders
decoders = assignTexNames(decoders);
decoders = assignColorNames(decoders);
decoders = assignColors(decoders);
decoders = assignPlotLineTypes(decoders);
decoders = assignBarColorStyles(decoders);

% and the kinematic variables
kinemat = assignKinTexNames(kinemat);
kinemat = assignLegendNames(kinemat);

% 
if exist('fracneurons','var'), varargout{1} = fracneurons; end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [B,p] = nearestRivalPlots(decoders,kinemat,binwidths,mkinds,monkeys,....
    standardname,rivalnames)
% R^2 and SNR: rEFH vs. nearest rival (\FigNearestRivalScatterPlots)

% Ns
Ndims = 2;                                  % two-dimensional space
Nstates = length(kinemat);
Norder = Nstates/Ndims;
Nbinwidths = length(binwidths);
Nmonkeys = length(mkinds);
Ndecoders = length(decoders);
Nrivals = length(rivalnames);


% "meta-data"
basefignum = 23;
monkeylabels    = repmat(monkeys,[Ndims,1]);
titleStrs       = reshape({kinemat.legendname},Ndims,[]);
texnamemat      = reshape({kinemat.texname},Ndims,[]);
filesuffices    = reshape({kinemat.name},Ndims,[]);

allRsqs = reshape(cat(4,decoders.CoD),[],Ndims,Norder,Nbinwidths,Ndecoders);
iStandard = strcmp({decoders(:).name},standardname);
rivals = find(ismember({decoders(:).name},rivalnames));


% for each binwidth
B = zeros(2,Nmonkeys,Norder,Nbinwidths,Nrivals);
p = zeros(1,Nmonkeys,Norder,Nbinwidths,Nrivals);
for iBinwidth = 1:Nbinwidths
    Nmsperbin = binwidths(iBinwidth);
            
    for iOrder = 1:Norder
        texnames = repmat(texnamemat(:,iOrder),[1,Nmonkeys]);        
        
        
        for iRival = 1:Nrivals
            iDecoder = rivals(iRival);
            
            scatterWithEqualityLine(...
                allRsqs(:,:,iOrder,iBinwidth,iDecoder),...
                allRsqs(:,:,iOrder,iBinwidth,iStandard),...
                mkinds,texnames,'CoD',0,1,...
                ['R$^2$, ',decoders(iDecoder).texname],...
                ['R$^2$, ',decoders(iStandard).texname],...
                titleStrs{1,iOrder},binwidths(iBinwidth),Nmsperbin,...
                decoders(iStandard).name,decoders(iDecoder).name,...
                monkeylabels,filesuffices{1,iOrder}(1:3),...
                basefignum+iOrder+10*iBinwidth+(iRival-1)*1000);
            
            [B(:,:,iOrder,iBinwidth,iRival), p(:,:,iOrder,iBinwidth,iRival)] =...
                scatterWithEqualityLine(...
                Rsq2SNR(allRsqs(:,:,iOrder,iBinwidth,iDecoder)),...
                Rsq2SNR(allRsqs(:,:,iOrder,iBinwidth,iStandard)),...
                mkinds,texnames,'SNR',0,10,...
                ['SNR (dB), ',decoders(iDecoder).texname],...
                ['SNR (dB), ',decoders(iStandard).texname],...
                titleStrs{1,iOrder},binwidths(iBinwidth),Nmsperbin,...
                decoders(iStandard).name,decoders(iDecoder).name,...
                monkeylabels,filesuffices{1,iOrder}(1:3),...
                basefignum+iOrder+10*iBinwidth+(iRival-1)*1000+100);
            
        end
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [beta1,p] = scatterWithEqualityLine(xx,yy,mkinds,varsymbols,...
    fitmetric,fitmetricMin,fitmetricMax,yrxlabel,yrylabel,varname,...
    binwidth,Nmsperbin,standardstr,rivalstr,monkeylabels,filesuffix,fignum)


% init
Nmonkeys = length(mkinds);
[~,machine] = system('hostname');
machine = strtrim(machine);
h = figure(fignum); clf; hold on;
clrs = [
    27,158,119;...      % greenish
    217,95,2;...        % orangish
    117,112,179;...     % purplish
    ]/255;
% cbrewer Dark2


% for each monkey...
beta1 = zeros(2,Nmonkeys);
p = zeros(1,Nmonkeys);
for iMk = 1:Nmonkeys
    
    xxx = xx(mkinds{iMk})';
    yyy = yy(mkinds{iMk})';
    Msessions = length(xxx);
    
    % H0: assume the slope is unity
    beta0 = [1; mean(yyy - xxx)];
    SS0 = sum( (yyy - [xxx,ones(Msessions,1)]*beta0).^2  );
    df0 = Msessions - 1; % intercept
    
    % H1: fit the slope
    [beta1(:,iMk),~,Yres] = linregress([xxx,ones(Msessions,1)],yyy);
    SS1 = sum(Yres.^2);
    df1 = Msessions - 2; % slope and intercept
    
    % f statistic
    f = ((SS0 - SS1)/(df0 - df1))/(SS1/df1);
    p(iMk) = 1 - fcdf(f,df0,df1);
    
    % scatter
    scatter(xx(mkinds{iMk},1),yy(mkinds{iMk},1),'x',...
        'MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
    scatter(xx(mkinds{iMk},2),yy(mkinds{iMk},2),'o',...
        'MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
    labels = arrayfun(@(ii)([monkeylabels{ii}, ', ',varsymbols{ii}]),...
        1:numel(monkeylabels),'UniformOutput',0);
end
plot([fitmetricMin,fitmetricMax],[fitmetricMin,fitmetricMax],'k');
hold off;

% annotate
legend(labels,'Location','SouthEast','Interpreter','Latex')
xlabel(yrxlabel)
ylabel(yrylabel)
title(sprintf('%s (%02dms bins)',varname,binwidth),'Interpreter','Latex');
axis equal tight

% if you're on a computer with LaTeX
switch machine
    case {'kobayashi-maru','CUPCAKE','Themistocles'}
        matlab2tikzWrapper(sprintf('SLMD/%s_%s_%0.3i_%s_vs_%s_scatter',...
            fitmetric,filesuffix,Nmsperbin,standardstr,rivalstr),h);
    otherwise
        fprintf('skipping tikz plots on this machine -- jgm\n');
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function SNR = Rsq2SNR(Rsq)

SNR = -10*log10(1 - Rsq);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function summaryStats(decoders,kinemat,binwidths,NdataTest,standardname,rivalnames)

% Ns
Ndims = 2;
Nstates = length(kinemat);
Norder = Nstates/Ndims;
Nsessions = size(NdataTest,1);
Nbinwidths = length(binwidths);

% indices
iStandard = strcmp({decoders(:).name},standardname);
rivals = find(ismember({decoders(:).name},rivalnames));


% for each "rival"
kinnames = reshape({kinemat.legendname},Ndims,[]);
for iRival = 1:length(rivals)
    iDecoder = rivals(iRival);
    
    
    keyboard
    
    % SNR
    [MeanImprovementSNR,ISSIGgauss,ISSIGperm] = signify(...
        Rsq2SNR(decoders(iStandard).CoD),Rsq2SNR(decoders(iDecoder).CoD),...
        NdataTest,0.01);
    printmatJGM(squeeze(MeanImprovementSNR)',...
        sprintf('Improvement (rEFH over %s) in SNR (dB)',...
        decoders(iDecoder).texname),...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',kinemat.name))
    printmatJGM(squeeze(ISSIGgauss)',...
        sprintf('Significant difference with %s? (rEFH=+1 over %s)',...
        decoders(iDecoder).texname),...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',kinemat.name))
    
    
    
    [MeanImprovementSNR,ISSIGgauss,ISSIGperm] = signify(...
        reshape(Rsq2SNR(decoders(iStandard).CoD),Nsessions*Ndims,Norder,Nbinwidths,[]),...
        reshape(Rsq2SNR(decoders(iDecoder).CoD),Nsessions*Ndims,Norder,Nbinwidths,[]),...
        repmat(NdataTest,[Ndims,1,1,1]),0.01);
    printmatJGM(squeeze(MeanImprovementSNR)',...
        sprintf('Improvement (rEFH over %s) in SNR (dB)',...
        decoders(iDecoder).texname),...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',kinnames{1,:}))
    printmatJGM(squeeze(ISSIGgauss)',...
        sprintf('Significant difference with %s? (rEFH=+1 over %s)',...
        decoders(iDecoder).texname),...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',kinnames{1,:}))
    
    
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [MeanImprovement,ISSIGgauss,ISSIGperm] =...
    signify(fitsStandard,fitsRival,NdataTest,alp)

% Ns
[Nsessions,Nstates,Nbins] = size(fitsStandard);
Nperms = 50000;


% based on the assumption of a normal distribution
DiffWithRival = fitsStandard - fitsRival;
MeanImprovement = sum(NdataTest.*DiffWithRival,1)./sum(NdataTest,1);
VrncImprovement = sum(NdataTest.*...
    (DiffWithRival - MeanImprovement).^2)./sum(NdataTest);
StdErrorOfMeanImprovementSNR = sqrt(VrncImprovement/Nsessions);
ISSIGgauss = (norminv(alp/2,...
    abs(MeanImprovement),StdErrorOfMeanImprovementSNR) > 0).*...
    sign(MeanImprovement);


% based on a permutation test
wtdFitsStandard = fitsStandard.*NdataTest./sum(NdataTest);
wtdFitsRival = fitsRival.*NdataTest./sum(NdataTest);

% cat rival first because you use diff below
bothFits = cat(1,wtdFitsRival,wtdFitsStandard);

% now sample
inds = ceil(2*Nsessions*rand(Nsessions*2*Nperms,1));
bootFits = reshape(bothFits(inds,:,:),Nsessions,2,Nperms,Nstates,Nbins);
bootDstrb = permute(diff(sum(bootFits),[],2),[3,4,5,1,2]);
pBoot = mean(MeanImprovement < bootDstrb);
ISSIGperm = pBoot < alp;


keyboard
bothFits = cat(1,fitsRival,fitsStandard);
MdataTest = repmat(NdataTest,[2,1,1]);
wts = permute(MdataTest./sum(MdataTest),[1,3,2]);
inds = categorsmpl(wts(:,1),Nsessions*2*Nperms,'IndexBased');
%%% might want to decouple across bins, although it's not necessary
bootFits = reshape(bothFits(inds,:,:),Nsessions,2,Nperms,Nstates,Nbins);
bootDstrb = permute(mean(diff(bootFits,[],2)),[3,4,5,1,2]);
pBoot = mean(MeanImprovement < bootDstrb);
ISSIGperm = pBoot < alp;


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function allDecodersBarPlots(decoders,kinemat,binwidths,NdataTest)
% All decoders summary bar plot

% all results together
allRsqs = cat(4,decoders.CoD);
iStatic = strcmp({decoders(:).name},'static');

% R^2
DiffWithStatic = allRsqs(:,:,:,~iStatic) - allRsqs(:,:,:,iStatic);
barPlotCore(DiffWithStatic,kinemat,binwidths,NdataTest,...
    'R$^2$ improvement',{decoders(~iStatic).barcolor},'CoD',...
    'extraTikzText','[baseline,trim axis left,trim axis right]');

% SNR
DiffWithStatic = Rsq2SNR(allRsqs(:,:,:,~iStatic)) -...
    Rsq2SNR(allRsqs(:,:,:,iStatic));
barPlotCore(DiffWithStatic,kinemat,binwidths,NdataTest,...
    'SNR (dB) improvement',{decoders(~iStatic).barcolor},'SNR',...
    'extraTikzText','[baseline,trim axis left,trim axis right]');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function barPlotCore(fitscores,kinemat,binwidths,NdataTest,...
    fitmetric,barcolors,fileprefix,varargin)

% Ns
[Nsessions,Nstates,Nbinwidths,Nfilters] = size(fitscores);

% plot the averages, with error bars
MeanFits = sum(NdataTest.*fitscores,1)./sum(NdataTest,1);
VrncFits = sum(NdataTest.*(fitscores - MeanFits).^2,1)./sum(NdataTest,1);
for iBinwidth = 1:Nbinwidths
    muDiff = permute(MeanFits(1,:,iBinwidth,:),[2,4,1,3]);
    stdErr = sqrt(VrncFits(1,:,iBinwidth,:)/Nsessions);
    stdErr = permute(stdErr,[2,1,4,3]);
    
    ymin = defaulter('ymin', min(0,min(vect(muDiff - stdErr))),varargin{:});
    tikzBarGraph(...
        1:Nstates,...
        muDiff,...
        cat(2,-stdErr,stdErr),...
        ymin,...
        {kinemat.texname},'kinematic variable',...
        sprintf('%s',fitmetric),'',...
        repmat(barcolors',[1,Nstates])',...
        2,0.32,'grouped',{},...
        sprintf('SLMD/%s_allvars_%03d_allmodels_bar',...
        fileprefix,binwidths(iBinwidth)),varargin{:});
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfBinwidthPlots(decoders,kinemat,binwidths,NdataTest,...
    SEMILOGPLOT)
% Within the .tex file, you may want to:
%   (1) bold the x-axis label 64 (for posters)


% Ns
Ndims = 2;                                  % two-dimensional space

% reshape
allRsqs = cat(4,decoders.CoD); 
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,binwidths,NdataTest,allRsqs,'CoD',0,0.8,...
    titleStrs,filesuffices,'binwidth','bin width (ms)','R$^2$',233,SEMILOGPLOT);

% SNR
effectOfCore(decoders,binwidths,NdataTest,Rsq2SNR(allRsqs),'SNR',0,6,...
    titleStrs,filesuffices,'binwidth','bin width (ms)','SNR (dB)',333,SEMILOGPLOT);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfNumberOfTrainingSamples(decoders,monkeys,SEMILOGPLOT)

% load data
[decoders,kinemat,~,trainingtimes,NdataTest] =...
    assembleData(decoders,monkeys,'Rsqs_064msBins_allTrainingTimes');

% reshape
Ndims = 2;
allRsqs = cat(4,decoders.CoD); 
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,trainingtimes,NdataTest,allRsqs,'CoD',0,0.8,...
    titleStrs,filesuffices,'trainingtime','seconds of training data',...
    'R$^2$',433,SEMILOGPLOT);

% SNR
effectOfCore(decoders,trainingtimes,NdataTest,Rsq2SNR(allRsqs),'SNR',0,6,...
    titleStrs,filesuffices,'trainingtime','seconds of training data',...
    'SNR (dB)',533,SEMILOGPLOT);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfNumberOfNeurons(decoders,kinemat,binwidths,Nneurons,...
    mkinds,monkeys,SEMILOGPLOT)


% useful function
centerByMonkey = @(XX)(cell2mat(arrayfun(@(ii)(...
    (XX(mkinds{ii},:,:,:) - mean(XX(mkinds{ii},:,:,:)))),...
    1:length(mkinds),'UniformOutput',0)'));

% first see if Nneurons correlated with performance in the canonical data
allRsqs = cat(4,decoders.CoD);
[Nsessions,Nstates,~,Nfilters] = size(allRsqs);

% you have to center the data before they're comparable across animals
NneuronsCtrd = centerByMonkey(Nneurons);
reconQualityCtrd = centerByMonkey(allRsqs);

yrcorrs = reshape(...
    corr(NneuronsCtrd,reshape(reconQualityCtrd(:,:,binwidths==64,:),...
    [Nsessions,Nstates*Nfilters])),[Nstates,Nfilters]);
printmatJGM(yrcorrs,...
    sprintf('Correlations between number of neurons and Rsq'),...
    sprintf('%s ',kinemat.name),sprintf('%s ',decoders.name))
fitmetric = 'R$^2$';
tikzBarGraph(...
    1:Nstates,...
    yrcorrs,...
    nan(Nstates,2,Nfilters,2),...
    min(0,min(vect(yrcorrs))),...
    {kinemat.texname},'kinematic variable',...
    sprintf('Corr(%s,N$_\\text{neurons}$)',fitmetric),'',...
    repmat({decoders.barcolor}',[1,Nstates])',...
    2,0.32,'grouped',{},...
    sprintf('SLMD/%s_allvars_%03d_allmodels_bar',...
    'corr_CoD_Nneurons',binwidths(3)));

reconQualityCtrd = centerByMonkey(Rsq2SNR(allRsqs));
yrcorrs = reshape(corr(...
    NneuronsCtrd,reshape(reconQualityCtrd(:,:,binwidths==64,:),...
    [Nsessions,Nstates*Nfilters])),[Nstates,Nfilters]);
printmatJGM(yrcorrs,...
    sprintf('Correlations between number of neurons and SNR'),...
    sprintf('%s ',kinemat.name),sprintf('%s ',decoders.name))
fitmetric = 'SNR';
tikzBarGraph(...
    1:Nstates,...
    yrcorrs,...
    nan(Nstates,2,Nfilters,2),...
    min(0,min(vect(yrcorrs))),...
    {kinemat.texname},'kinematic variable',...
    sprintf('Corr(%s,N$_\\text{neurons}$)',fitmetric),'',...
    repmat({decoders.barcolor}',[1,Nstates])',...
    2,0.32,'grouped',{},...
    sprintf('SLMD/%s_allvars_%03d_allmodels_bar',...
    'corr_SNR_Nneurons',binwidths(3)));



% load data
[decoders,kinemat,binwidths,trainingtimes,NdataTest,~,~,fracneurons] =...
    assembleData(decoders,monkeys,'Rsqs_064msBins_fewerneurons');

% reshape
Ndims = 2;
allRsqs = cat(4,decoders.CoD);
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,fracneurons,NdataTest,allRsqs,'CoD',0,0.8,...
    titleStrs,filesuffices,'numneurons','fraction of all neurons',...
    'R$^2$',433,SEMILOGPLOT);

% SNR
effectOfCore(decoders,fracneurons,NdataTest,Rsq2SNR(allRsqs),'SNR',0,6,...
    titleStrs,filesuffices,'numneurons','fraction of all neurons',...
    'SNR (dB)',533,SEMILOGPLOT);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfCore(decoders,xdata,NdataTest,fitscores,...
  fitmetric,fitmetricMin,fitmetricMax,titleStrs,filesuffices,sweptparam,...
  yrxlabel,yrylabel,basefignum,SEMILOGPLOT)


% Ns
Ndecoders = length(decoders);
Nx = length(xdata);
[Ndims,Norder] = size(titleStrs);

% for plotting
[~,machine] = system('hostname');
machine = strtrim(machine);
clrNameClrPairs = arrayfun(@(ii)({decoders(ii).colorname,decoders(ii).color}),...
    1:Ndecoders,'UniformOutput',false);
if SEMILOGPLOT, xdataplot = log2(xdata); else, xdataplot = xdata; end

% average and reshape
MeanFits = permute(...
    mean(reshape(sum(NdataTest.*fitscores,1)./sum(NdataTest,1),...
    [],Ndims,Norder,Nx,Ndecoders),2),[3,4,5,1,2]);

% loop through pos, vel, acc
for iOrder = 1:Norder
    filesuffix = filesuffices{1,iOrder}(1:3);
    h = figure(basefignum + iOrder); clf; hold on;
    for iDecoder = 1:Ndecoders
        plot(xdataplot,MeanFits(iOrder,:,iDecoder),...
            decoders(iDecoder).linetype,...
            'Color',decoders(iDecoder).color,...
            'linewidth',2.0);
    end
    set(gca,'Xtick',xdataplot);
    if SEMILOGPLOT
        axis([xdataplot(1),xdataplot(end),fitmetricMin,fitmetricMax]);
    else
        axis([0,xdataplot(end),fitmetricMin,fitmetricMax]);
    end
    set(gca,'Xticklabels',xdata);
    
    hold off;
    title(titleStrs{1,iOrder})
    legend({decoders.texname},'Location','SouthEast')
    xlabel(yrxlabel);
    ylabel(yrylabel);
    switch machine
        case {'kobayashi-maru','CUPCAKE','Themistocles'}
            matlab2tikzWrapper(...
                sprintf('SLMD/%s_%s_vs_%s_allmodels_line',...
                fitmetric,filesuffix,sweptparam),h,...
                'extraColors',clrNameClrPairs,...
                'extraTikzpictureOptions',...
                {'baseline','trim axis left','trim axis right'});
        otherwise
            fprintf('skipping tikz plots on this machine -- jgm\n');
    end
    
    
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfREFHcomponents(monkeys)

% % the static decoder--a baseline
% decoderStatic.name = 'static';
% [decoderStatic,kinematPois,binwidthsPois] =...
%     assembleData(decoderStatic,monkeys,'Rsqs_allBins');
decoderUKF.name = 'ukf';
decoderUKF = assembleData(decoderUKF,monkeys,'Rsqs_allBins');


% re-collect the REFH (Poisson) decoders
[decodersPois(1:2).name] = deal('refhstatic','refhdynamic');
[decodersPois,kinematPois,binwidthsPois,~,NdataTestPois] =...
    assembleData(decodersPois,monkeys,'Rsqs_allBins');

% the REFH, StandardNormal decoders
[decodersStdNrml(1:2).name] = deal('refhstatic_stdnrml','refhdynamic_stdnrml');
[decodersStdNrml,kinematStdNrml,binwidthsStdNrml,~,NdataTestStdNrml] =...
    assembleData(decodersStdNrml,monkeys,'Rsqs_allBins_StandardNormal');

% the rEFH, decoding only from hidden units, and regularization
[decodersHidsOnly(1:2).name] = deal('refhstatic_hidsonly','refhdynamic_hidsonly');
[decodersHidsOnly,kinematHidsOnly,binwidthsHidsOnly,~,NdataTestHidsOnly] =...
    assembleData(decodersHidsOnly,monkeys,'Rsqs_allBins_HidsOnly');


% now make sure these data are consistent with the Poisson data
kinemat = kinematPois;
if ~all(strcmp({kinematStdNrml.name},{kinemat.name}))
    error('mismatching kinematic variables for StandardNormal results')
end
if ~all(strcmp({kinematHidsOnly.name},{kinematPois.name}))
    error('mismatching kinematic variables for HidsOnly results')
end
binwidths = binwidthsPois;
if ~all(binwidthsStdNrml == binwidths)
    error('mismatching binwidths for StandardNormal results')
end
if ~all(binwidthsHidsOnly == binwidthsPois)
    error('mismatching binwidths for HidsOnly results')
end
NdataTest = NdataTestPois;
if ~all(vect(NdataTestStdNrml == NdataTest))
    error('mismatching number of testing data for StandardNormal results')
end
if ~all(vect(NdataTestHidsOnly == NdataTest))
    error('mismatching number of testing data for HidsOnly results')
end



% now plot
[decodersHidsOnly(:).barcolor] = deal('black!20!white');
[decodersStdNrml(:).barcolor] = deal('black!50!white');
[decodersPois(:).barcolor] = deal('black');
decoders = [decodersHidsOnly,decodersStdNrml,decodersPois];
staticdecoders = {'refhstatic_hidsonly','refhstatic_stdnrml','refhstatic'};
for iDecoder = 1:length(staticdecoders)
    staticdecoder = staticdecoders{iDecoder};
    decoders = cleverAssign(decoders,staticdecoder,'barcolor',...
        [decoders(strcmp({decoders(:).name},staticdecoder)).barcolor,...
        ',postaction={pattern=north east lines,pattern color=white}']);
end
allRsqs = cat(4,decoders.CoD);


% combine across dimensions
Ndims = 2; % fact 
tmp = reshape({kinemat.legendname},Ndims,[]);
[xaxislabels(1:length(kinemat)/Ndims).texname] = deal(tmp{1,:});

% R^2
DiffWithUKF = allRsqs - decoderUKF.CoD;
DiffWithUKF = reshape(DiffWithUKF,size(allRsqs,1)*Ndims,...
    length(kinemat)/Ndims,length(binwidths),[]);
barPlotCore(DiffWithUKF,xaxislabels,binwidths,...
    repmat(NdataTest,[Ndims,1,1]),'R$^2_\text{rEFH}$ - R$^2_\text{UKF}$',...
    {decoders.barcolor},'CoD_rEFH');

% SNR
bruteforcetext = getbruteforcetext;
DiffWithUKF = Rsq2SNR(allRsqs) - Rsq2SNR(decoderUKF.CoD);
DiffWithUKF = reshape(DiffWithUKF,size(allRsqs,1)*Ndims,...
   length(kinemat)/Ndims,length(binwidths),[]);
barPlotCore(DiffWithUKF,xaxislabels,binwidths,...
    repmat(NdataTest,[Ndims,1,1]),'SNR$_\text{rEFH}$ - SNR$_\text{UKF}$ (dB)',...
    {decoders.barcolor},'SNR_rEFH','extraAxisText',bruteforcetext,...
    'ymin',-0.5);





end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function kinemat = assignKinTexNames(kinemat)
% TeX names for kinematic variables

kinemat = cleverAssign(kinemat,'posx','texname','$x$');
kinemat = cleverAssign(kinemat,'posy','texname','$y$');
kinemat = cleverAssign(kinemat,'velx','texname','$\dot x$');
kinemat = cleverAssign(kinemat,'vely','texname','$\dot y$');
kinemat = cleverAssign(kinemat,'accx','texname','$\ddot x$');
kinemat = cleverAssign(kinemat,'accy','texname','$\ddot y$');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function kinemat = assignLegendNames(kinemat)

kinemat = cleverAssign(kinemat,'posx','legendname','Position');
kinemat = cleverAssign(kinemat,'posy','legendname','Position');
kinemat = cleverAssign(kinemat,'velx','legendname','Velocity');
kinemat = cleverAssign(kinemat,'vely','legendname','Velocity');
kinemat = cleverAssign(kinemat,'accx','legendname','Acceleration');
kinemat = cleverAssign(kinemat,'accy','legendname','Acceleration');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignTexNames(decoders)
% legend (order of these statements doesn't matter)

% standard decoders
decoders = cleverAssign(decoders,'static',      'texname','lin. regr');
decoders = cleverAssign(decoders,'kfobs',       'texname','KF$_\text{Obs}$');
decoders = cleverAssign(decoders,'kfemstatic',  'texname','KF$_\text{EM}$ (static)');
decoders = cleverAssign(decoders,'kfemdynamic', 'texname','KF$_\text{EM}$ (+KF)');
decoders = cleverAssign(decoders,'ukf',         'texname','UKF');
decoders = cleverAssign(decoders,'refhstatic',  'texname','rEFH (static)');
decoders = cleverAssign(decoders,'refhdynamic', 'texname','rEFH (+KF)');

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml',  'texname',...
    'rEFH$_\text{Gauss}$ (static)');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml', 'texname',...
    'rEFH$_\text{Gauss}$ (+KF)');
decoders = cleverAssign(decoders,'refhstatic_hidsonly',  'texname',...
    'rEFH$_\text{no reg.}$ (static)');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly', 'texname',...
    'rEFH$_\text{no reg.}$ (+KF)');


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignColorNames(decoders)
% color names (order of these statements doesn't matter)
%
% NB that this is actually more important than assignColors, since the
% color name written into the tex file by this macro will ultimately be
% defined externally in a master tex file, overriding the color definitions
% written on the basis of asignColor.

% standard decoders
decoders = cleverAssign(decoders,'static',      'colorname','xtraclr');
decoders = cleverAssign(decoders,'kfobs',       'colorname','OBSclr');
decoders = cleverAssign(decoders,'kfemstatic',  'colorname','EMclr');
decoders = cleverAssign(decoders,'kfemdynamic', 'colorname','EMclr');
decoders = cleverAssign(decoders,'ukf',         'colorname','pinkish');
decoders = cleverAssign(decoders,'refhstatic',  'colorname','EFHclr');
decoders = cleverAssign(decoders,'refhdynamic', 'colorname','EFHclr');

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml','colorname','EFHclr');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml','colorname','EFHclr');
decoders = cleverAssign(decoders,'refhstatic_hidsonly','colorname','gray');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','colorname','gray');


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignColors(decoders)
% color names (order of these statements doesn't matter)
%
% But see assignColorNames, above.

setColors;

% standard decoders
decoders = cleverAssign(decoders,'static',      'color',XTRAcolor);
decoders = cleverAssign(decoders,'kfobs',       'color',OBScolor);
decoders = cleverAssign(decoders,'kfemstatic',  'color',EMcolor);
decoders = cleverAssign(decoders,'kfemdynamic', 'color',EMcolor);
decoders = cleverAssign(decoders,'ukf',         'color',PINKish);
decoders = cleverAssign(decoders,'refhstatic',  'color',EFHcolor);
decoders = cleverAssign(decoders,'refhdynamic', 'color',EFHcolor);

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml','color',EFHcolor);
decoders = cleverAssign(decoders,'refhdynamic_stdnrml','color',EFHcolor);
decoders = cleverAssign(decoders,'refhstatic_hidsonly','color',[0.2,0.2,0.2]);
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','color',[0.2,0.2,0.2]);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignPlotLineTypes(decoders)
% the plot line types (order of these statements doesn't matter)

% standard decoders
decoders = cleverAssign(decoders,'static',      'linetype','-');
decoders = cleverAssign(decoders,'kfobs',       'linetype','-');
decoders = cleverAssign(decoders,'kfemstatic',  'linetype','--');
decoders = cleverAssign(decoders,'kfemdynamic', 'linetype','-');
decoders = cleverAssign(decoders,'ukf',         'linetype','-');
decoders = cleverAssign(decoders,'refhstatic',  'linetype','--');
decoders = cleverAssign(decoders,'refhdynamic', 'linetype','-');

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml','linetype','--');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml','linetype','-');
decoders = cleverAssign(decoders,'refhstatic_hidsonly','linetype','--');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','linetype','-');
%%% you would perhaps like to alternate colors....



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignBarColorStyles(decoders)
% the bar color styles (order of these statements doesn't matter)
% 
% The "static" decoders are hatched.

[decoders(:).barcolor] = deal(decoders(:).colorname);

%%% you *could* do this based on the 'linetype'---assuming
%%% assignPlotLineTypes has already been run.  Ideally you'd have an if
%%% statement and try both.....
staticdecoders = {'kfemstatic','refhstatic','refhstatic_stdnrml'};
%%%

for iDecoder = 1:length(staticdecoders)
    staticdecoder = staticdecoders{iDecoder};
    decoders = cleverAssign(decoders,staticdecoder,'barcolor',...
        [decoders(strcmp({decoders(:).name},staticdecoder)).barcolor,...
        ',postaction={pattern=north east lines,pattern color=white}']);
end

gaussiandecoders = {'refhdynamic_stdnrml','refhstatic_stdnrml'};
for iDecoder = 1:length(gaussiandecoders)
    gaussiandecoder = gaussiandecoders{iDecoder};
    decoders = cleverAssign(decoders,gaussiandecoder,'barcolor',...
        [decoders(strcmp({decoders(:).name},gaussiandecoder)).barcolor,...
        ',draw=black']);
end 

hidsonlydecoders = {'refhdynamic_hidsonly','refhstatic_hidsonly'};
for iDecoder = 1:length(hidsonlydecoders)
    hidsonlydecoder = hidsonlydecoders{iDecoder};
    decoders = cleverAssign(decoders,hidsonlydecoder,'barcolor',...
        [decoders(strcmp({decoders(:).name},hidsonlydecoder)).barcolor,...
        ',draw=black']);
end 




end
%-------------------------------------------------------------------------%

% %-------------------------------------------------------------------------%
% function decoders = assignCoDs(RsqStatic,RsqKF_Obs,RsqKF_EMstatic,...
%     RsqKF_EMdynamic,RsqUKF,RsqREFHstatic,RsqREFHdynamic,decoders)
% 
% % standard decoders
% decoders = cleverAssign(decoders,'static',      'CoD',RsqStatic);
% decoders = cleverAssign(decoders,'kfobs',       'CoD',RsqKF_Obs);
% decoders = cleverAssign(decoders,'kfemstatic',  'CoD',RsqKF_EMstatic);
% decoders = cleverAssign(decoders,'kfemdynamic', 'CoD',RsqKF_EMdynamic);
% decoders = cleverAssign(decoders,'ukf',         'CoD',RsqUKF);
% decoders = cleverAssign(decoders,'refhstatic',  'CoD',RsqREFHstatic);
% decoders = cleverAssign(decoders,'refhdynamic', 'CoD',RsqREFHdynamic);
% 
% % non-standard decoders
% decoders = cleverAssign(decoders,'refhstatic_stdnrml','CoD',RsqREFHstatic);
% decoders = cleverAssign(decoders,'refhdynamic_stdnrml','CoD',RsqREFHdynamic);
% %%% Notice that these won't be invoked unless refhstatic_stdnrml is an
% %%% element in decoders
% 
% end
% %-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function structarray = cleverAssign(structarray,elementname,fieldname,fieldcontent)
% This clever little one-liner finds the element elementname (e.g., 'posx')
% in the structure array structarray (e.g., kinemat), and assigns the value
% of fieldcontent (e.g., 'Position') to the field fieldname (e.g. texname).
%
% Each element of the array must therefore have a field called 'name' in 
% which the elementname is set.  But what makes this fxn clever is that it 
% doesn't fail when the element name does not exist.

[structarray(strcmp({structarray(:).name},elementname)).(fieldname)] =...
    deal(fieldcontent);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function bruteforcetext = getbruteforcetext

bruteforcetext = sprintf([...
    'restrict y to domain*=\\pgfkeysvalueof{/pgfplots/ymin}:',...
    '\\pgfkeysvalueof{/pgfplots/ymax}, %% Cut values off at ymin\n',...
    'visualization depends on=rawy\\as\\rawy, %% Save the unclipped values\n',...
    'after end axis/.code={ %% Draw line indicating break\n',...
    '\t\\draw [ultra thick, white, decoration={snake, amplitude=1pt}, decorate] ',...
    '(rel axis cs:0,0.1) -- (rel axis cs:1,0.1);\n',...
    '},\n',...
    'clip=false,\n',...
    'nodes near coords={%%\n',...
    '\t\\pgfkeys{/pgf/fpu=true} %% Switch on the fpu library\n',...
    '\t\\pgfmathparse{\\rawy<\\pgfkeysvalueof{/pgfplots/ymin}} %% Do the comparison\n',...
    '\t\\pgfmathfloattofixed{\\pgfmathresult} %% convert the result to fixed point\n',...
    '\t\\pgfkeys{/pgf/fpu=false} %% switch off the fpu library\n',...
    '\t\\ifdim\\pgfmathresult pt=1pt %% If the condition was true...\n',...
    '\t\t\\pgfmathprintnumber{\\rawy}\n',...
    '\t\\else\n',...
    '\t\t%%\n',...
    '\t\\fi\n',...
    '},\n',...
    'every node near coord/.append style={yshift =-0.1in,font=\\footnotesize},\n',...
    'every x tick label/.append style={yshift = -0.2in},\n',...
    ]);


%every node near coord/.append style={rotate=90, anchor=east, xshift =-0.15in,font=\footnotesize},

end
%-------------------------------------------------------------------------%


