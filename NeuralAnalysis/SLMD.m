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
% Revised: 09/30/17
%   -added plotIllustrativeReconstructions
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
[decoders(1:8).name] = deal('static','kfobs','kfemstatic',...
    'kfemdynamic','wf','ukf','refhstatic','refhdynamic');

% load the saved results
[decoders,kinemat,binwidths,trainingtimes,NdataTest,Nneurons,mkinds] =...
    assembleData(decoders,monkeys,'Rsqs_BinwidthSweep_Poisson');
fprintf('"canonical" training time is %d seconds\n',trainingtimes);


% you don't usually want to plot everything....
keyboard



% plot an illustrative decoding example (SNRs for refhdynamic for session
%   21 are close to their avg.)
plot_binwidth = 64;
plot_matrix = {'t', 'posy';
    't', 'vely';
    't', 'accy'};
decoders_with_lags = plotIllustrativeReconstructions(plot_matrix,'Indy',...
    21,plot_binwidth,320,1.0,'binwidths','Poisson',kinemat,...
    {'groundtruth','kfobs','kfemdynamic','ukf','refhdynamic'},...
    [60, 75]);
plotReconstructionLaggedSpikeCoDs(decoders_with_lags,kinemat,plot_binwidth)


% rEFH vs. nearest rival
[BByMk,pByMk,BAllMks,pAllMks] = nearestRivalPlots(decoders,kinemat,...
    binwidths,mkinds,monkeys,'refhdynamic',{'kfobs','ukf'});

% how much improvement for pos, vel, acc?
summaryStats(decoders,kinemat,binwidths,NdataTest,'refhdynamic',...
    {'kfobs','ukf'},0.01,'ms');

% bar summary plots (\FigAllDecodersBarPlots)
allDecodersBarPlots(decoders,kinemat,binwidths,NdataTest,'alldecoders');

% effect of binwidths
effectOfBinwidthPlots(decoders,kinemat,binwidths,NdataTest,0,...
    'alldecoders',[0,0.8],[0,6]);

% effect of number of training samples
effectOfNumberOfTrainingSamples(decoders,monkeys,0,...
    'alldecoders',[0,0.8],[0,6]);

% effect of dropping out neurons
effectOfNumberOfNeurons(decoders,kinemat,binwidths,NdataTest,Nneurons,...
    mkinds,monkeys,0,'alldecoders',[0,0.8],[0,6]);

% effect of rEFH components
effectOfREFHcomponents(monkeys);

% using tapped decoders
tappedDecoding(monkeys)




% for the replies to the reviewers/supplemental
keyboard
compare_LGDS_with_Fixed_20_Hids;
compare_EM4LDS_jgm_vs_zg;
allDecodersBarPlots_one_monkey('Indy')
allDecodersBarPlots_one_monkey('Loco')
effectOfNumberOfNeurons(decoders,kinemat,binwidths,NdataTest,Nneurons,...
    mkinds,monkeys,0,'alldecoders_negative_cods',[0,0.8],[-2,6]);
effectOfNumberOfNeurons_scatter(decoders,kinemat,binwidths,...
    Nneurons,mkinds,monkeys,{'refhstatic','refhdynamic','kfemstatic','ukf'});
effectsOfSparsityTargetandNumberOfHiddenUnits(monkeys);
%-------------------------------------------------------------------------%



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
    load(sprintf('%sRBMish/BMI/%s_%s.mat',getdir('data'),fileprefix,monkeys{iMk}));
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
    try
        [~,iDecoders] = ismember({decoders(:).name},decodernames);
        %%%allRsqs = cat(1,allRsqs,Rsqs(:,:,:,iDecoders));
        allRsqs = cat(1,allRsqs,Rsqs(:,:,:,iDecoders,:));
    catch ME
        keyboard
    end
    
    
    % accumulate other useful numbers
    Msessions       = size(Rsqs,1);
    mkinds{iMk}     = Nsessions + (1:Msessions);
    Nsessions       = Nsessions + Msessions;
    NdataTest_all   = [NdataTest_all; NdataTest];
    if exist('Nneurons','var'), Nneurons_all=[Nneurons_all;Nneurons]; end
    
    
end

% now break apart the big tensor to assign CoDs to each decoder
for iDecoder = 1:length(decoders)
    %%%decoders(iDecoder).CoD = allRsqs(:,:,:,iDecoder);
    decoders(iDecoder).CoD = allRsqs(:,:,:,iDecoder,:);
end

% also assign the other useful bits for these decoders
decoders = assignTexNames(decoders);
decoders = assignColorNames(decoders);
decoders = assignColors(decoders);
decoders = assignPlotLineTypes(decoders);
decoders = assignBarColorStyles(decoders);
decoders = assignMarks(decoders);

% and the kinematic variables
kinemat = assignKinTexNames(kinemat);
kinemat = assignLegendNames(kinemat);

%
if exist('fracneurons','var'), varargout{1} = fracneurons; end
if exist('hidsensoryratios','var'), varargout{2} = hidsensoryratios; end
if exist('phidtargets','var'), varargout{3} = phidtargets; end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [BByMk,pByMk,BAllMks,pAllMks] = nearestRivalPlots(...
    decoders,kinemat,binwidths,mkinds,monkeys,standardname,rivalnames)
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
BByMk = zeros(2,Nmonkeys,Norder,Nbinwidths,Nrivals);
pByMk = zeros(1,Nmonkeys,Norder,Nbinwidths,Nrivals);
BAllMks = zeros(2,Norder,Nbinwidths,Nrivals);
pAllMks = zeros(1,Norder,Nbinwidths,Nrivals);

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
            
            [BByMk(:,:,iOrder,iBinwidth,iRival), pByMk(:,:,iOrder,iBinwidth,iRival),...
                BAllMks(:,iOrder,iBinwidth,iRival), pAllMks(:,iOrder,iBinwidth,iRival)] = ...
                scatterWithEqualityLine(...
                Rsq2SNR(allRsqs(:,:,iOrder,iBinwidth,iDecoder)),...
                Rsq2SNR(allRsqs(:,:,iOrder,iBinwidth,iStandard)),...
                mkinds,texnames,'SNR',-1,10,...
                ['SNR (dB), ',decoders(iDecoder).texname],...
                ['SNR (dB), ',decoders(iStandard).texname],...
                titleStrs{1,iOrder},binwidths(iBinwidth),Nmsperbin,...
                decoders(iStandard).name,decoders(iDecoder).name,...
                monkeylabels,filesuffices{1,iOrder}(1:3),...
                basefignum+iOrder+10*iBinwidth+(iRival-1)*1000+100);
            
        end
    end
end


% check if slopes are significantly different from 1
for iRival = 1:Nrivals
    iDecoder = rivals(iRival);
    
    % the 1 means slope rather than intercept
    printmatJGM(squeeze(BAllMks(1,:,:,iRival))',...
        sprintf('SNR regression slopes, %s vs. %s',...
        decoders(iStandard).texname,decoders(iDecoder).texname),...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',titleStrs{1,:}))
    
    printmatJGM(squeeze(pAllMks(:,:,:,iRival))',...
        'p value for being different from unity-slope model',...
        sprintf('%dms ',binwidths),...
        sprintf('%s ',titleStrs{1,:}))
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [BByMk,pByMk,BAllMks,pAllMks] =...
    scatterWithEqualityLine(xx,yy,mkinds,varsymbols,...
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
BByMk = zeros(2,Nmonkeys);
pByMk = zeros(1,Nmonkeys);
for iMk = 1:Nmonkeys
    
    % is the advantage improving with SNR?  cf. w/model with m=1 (each mk)
    xxx = vect(xx(mkinds{iMk},:));
    yyy = vect(yy(mkinds{iMk},:));
    Msessions = length(xxx);
    
    % H0: assume the slope is unity
    beta0 = [1; mean(yyy - xxx)];
    SS0 = sum( (yyy - [xxx,ones(Msessions,1)]*beta0).^2  );
    df0 = Msessions - 1; % intercept
    
    % H1: fit the slope
    [BByMk(:,iMk),~,Yres] = linregress([xxx,ones(Msessions,1)],yyy);
    SS1 = sum(Yres.^2);
    df1 = Msessions - 2; % slope and intercept
    
    % f statistic
    f = ((SS0 - SS1)/(df0 - df1))/(SS1/df1);
    pByMk(iMk) = 1 - fcdf(f,df0,df1);
    % scatter
    scatter(xx(mkinds{iMk},1),yy(mkinds{iMk},1),150,'^',...
        'MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
    scatter(xx(mkinds{iMk},2),yy(mkinds{iMk},2),'o',...
        'MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
    labels = arrayfun(@(ii)([monkeylabels{ii}, ', ',varsymbols{ii}]),...
        1:numel(monkeylabels),'UniformOutput',0);
    
end
plot([fitmetricMin,fitmetricMax],[fitmetricMin,fitmetricMax],'k');

% is the advantage improving with SNR?  cf. w/model with m=1 (all monkeys)
yyy = yy(:);
xxx = xx(:);
Msessions = length(xxx);

% H0: assume the slope is unity
beta0 = [1; mean(yyy - xxx)];
SS0 = sum( (yyy - [xxx,ones(Msessions,1)]*beta0).^2  );
df0 = Msessions - 1; % intercept

% H1: fit the slope
[BAllMks,~,Yres] = linregress([xxx,ones(Msessions,1)],yyy);
SS1 = sum(Yres.^2);
df1 = Msessions - 2; % slope and intercept
plot([fitmetricMin,fitmetricMax],[fitmetricMin,fitmetricMax]*BAllMks(1) + BAllMks(2),'k--');

hold off;

% f statistic
f = ((SS0 - SS1)/(df0 - df1))/(SS1/df1);
pAllMks = 1 - fcdf(f,df0,df1);




% annotate
legend(labels,'Location','SouthEast','Interpreter','Latex')
xlabel(yrxlabel)
ylabel(yrylabel)
titlestr = sprintf('%s (%02dms bins)',varname,binwidth);
title(titlestr,'Interpreter','Latex');
%ax = axis;
%axis([ax(1), ax(2), fitmetricMin,fitmetricMax])
axis equal % tight
axis([fitmetricMin,fitmetricMax,fitmetricMin,fitmetricMax])


% say which points have been cut off the plot
offplot = find(...
    (yyy > fitmetricMax)|(yyy < fitmetricMin)|...
    (xxx > fitmetricMax)|(xxx < fitmetricMin));
for iExcluded = 1:length(offplot)
    fprintf('%s (%s): not plotting point at (%.3g,%.3g)\n',...
        titlestr,fitmetric,xxx(offplot(iExcluded)),yyy(offplot(iExcluded)));
end
fprintf('\n');


% if you're on a computer with LaTeX
switch machine
    case {'kobayashi-maru','CUPCAKE','Themistocles'}
        clrNameClrPairs = arrayfun(@(ii)(...
            {[lower(monkeylabels{1,ii}),'clr'],clrs(ii,:)}),...
            1:size(clrs,2),'UniformOutput',false);
        matlab2tikzWrapper(sprintf('SLMD/%s_%s_%0.3i_%s_vs_%s_scatter',...
            fitmetric,filesuffix,Nmsperbin,standardstr,rivalstr),h,...
            'extraColors',clrNameClrPairs);
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
function decoders = plotIllustrativeReconstructions(plot_matrix,...
    monkey,iSession,binwidth,trainingtime,fracneuron,swept_param,...
    obsv_dstrb,kinemat,decoder_names,plot_range_seconds)

path(path,[getdir('code'),'UKF2_code_pack']);

% Ns
Nmspers         = 1000; % fact
base_fignum     = 100;
[Nplots,Naxes]  = size(plot_matrix);

% other useful things
t = plot_range_seconds(1):(binwidth/Nmspers):plot_range_seconds(2);
bins = ceil(t./(binwidth/1000));

% pull out the decoders you care about
% decoder_inds = arrayfun( @(iDecoderName)(find(strcmp({all_decoders(:).name},...
%     decoder_names{iDecoderName}))),1:length(decoder_names));

[yr_decoders(1:length(decoder_names)).name] = deal(decoder_names{:});
yr_decoders = assignTexNames(yr_decoders);
yr_decoders = assignColorNames(yr_decoders);
yr_decoders = assignColors(yr_decoders);
yr_decoders = assignPlotLineTypes(yr_decoders);
yr_decoders = assignMarks(yr_decoders);

[~,params,refh_full_path] = loadEFHspikecounts(monkey,iSession,...
    binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,'EFH');

% convert matrix of plot strings to indices
decoder_names = {kinemat.name, 't'}; % include t! a kind of hack
plot_matrix_flattened = plot_matrix(:);
kinemat_plot_inds = reshape(bsxfun(@(ii,jj)(...
    strcmp( plot_matrix_flattened(ii), decoder_names(jj) )),...
    (1:length(plot_matrix_flattened))', 1:length(decoder_names) )*...
    (1:length(decoder_names))', [Nplots,Naxes]);

% init plotter
for iPlot = 1:Nplots, h(iPlot) = figure(base_fignum+iPlot); clf; hold on; end

% plot reconstructions
[SNR,decoders] = plot_reconstructions_core(yr_decoders,...
    monkey,kinemat_plot_inds,params,t,bins,length(kinemat),iSession,...
    binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,...
    refh_full_path,base_fignum);

% overlay plot of actual movements, label/legend/etc.
figure_name_prefix = sprintf(['%s_%s_%s_%imsBins_%isTrainingtime_',...
    '%iFracneurons'],'SLMD/illustrativeReconstructions',obsv_dstrb,...
    monkey,binwidth,trainingtime,fracneuron);
finalize_reconstruction_plots(SNR,kinemat_plot_inds,kinemat,yr_decoders,...
    figure_name_prefix,h,base_fignum);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [SNR,yr_decoders] = plot_reconstructions_core(yr_decoders,...
    monkey,kinemat_plot_inds,params,t,bins,Nstates,iSession,...
    binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,...
    refh_complete_path,base_fignum)

% Ns
Ndecoders   = length(yr_decoders);
Nplots      = size(kinemat_plot_inds,1);
Nslag       = 3;

% load data
if checkGPUavailability, dataclass='gpuArray'; else, dataclass='double';end
[Xtrain,Qtrain] = params.getLatents([],dataclass,...
    'sequencelength','singlesequence');
[Rtrain,Qtrain] = params.getData(Xtrain,Qtrain);
[Rtest,Xtest,Qtest] = params.getTestData(dataclass);
SStot = sum((Xtest - mean(Xtest)).^2);


% now get decoded trajs and plot them
SNR = zeros(Ndecoders,Nstates+1); % +1 for time
for iDecoder = 1:Ndecoders
    switch yr_decoders(iDecoder).name
        case 'groundtruth'
            this_Xhat = Xtest;
            this_Rsqs = [];
        case {'refhstatic','refhdynamic'}
            if ~exist('wts','var')
                load(refh_complete_path,'wts');
                [~,Xhat_refhstatic,Rsq_refhstatic,...
                    Xhat_refhdynamic,Rsq_refhdynamic] =...
                    linearREFHDecoder(Rtrain,Xtrain,Qtrain,...
                    Rtest,Xtest,Qtest,wts,params,SStot);
                clear wts
                
                %%% hack: you want this in any case
                if ~any(ismember({yr_decoders.name},'refhstatic'))
                    yr_decoders(end+1).lagCoD =...
                        getLagRsqs(Xhat_refhstatic,Rtest,Nslag,binwidth);
                    yr_decoders(end).name = 'refhstatic';
                    yr_decoders(end) = assignTexNames(yr_decoders(end));
                    yr_decoders(end) = assignColorNames(yr_decoders(end));
                    yr_decoders(end) = assignColors(yr_decoders(end));
                    yr_decoders(end) = assignPlotLineTypes(yr_decoders(end));
                    yr_decoders(end) = assignMarks(yr_decoders(end));
                end
                %%%
            end
            if strcmp(yr_decoders(iDecoder).name, 'refhstatic')
                this_Xhat = Xhat_refhstatic; clear Xhat_refhstatic
                this_Rsqs = Rsq_refhstatic;
            else
                this_Xhat = Xhat_refhdynamic; clear Xhat_refhdynamic
                this_Rsqs = Rsq_refhdynamic;
            end
            
        case {'kfemstatic','kfemdynamic'}
            
            if ~exist('LDSparamsEM','var')
                LDSparamsEM = loadEFHspikecounts(monkey,iSession,...
                    binwidth,trainingtime,fracneuron,swept_param,...
                    obsv_dstrb,'LDS');                
                
                % train and test
                [LDSparamsEM,~,~,Xhat_kfemstatic,Rsq_kfemstatic,...
                    Xhat_kfemdynamic,Rsq_kfemdynamic] =....
                    latentStateLDSDecoder([],Rtrain,Xtrain,...
                    Rtest,Xtest,params.Nmsperbin,[],LDSparamsEM,SStot);
                
            end
            if strcmp(yr_decoders(iDecoder).name,'kfemstatic')
                this_Xhat = Xhat_kfemstatic; clear Xhat_kfemstatic
                this_Rsqs = Rsq_kfemstatic;
            else
                this_Xhat = Xhat_kfemdynamic; clear Xhat_kfemdynamic
                this_Rsqs = Rsq_kfemdynamic;
            end
            
        case 'ukf'
            %%% would be better to load something saved....
            [~,this_Xhat,this_Rsqs] = UKFDecoder(Rtrain,Xtrain,Qtrain.T,...
                Rtest,Xtest,SStot,'ftaps',0,'ptaps',1,'htaps',0);
        case 'kfobs'
            [~,this_Xhat,this_Rsqs] = fullyObservedLDSDecoder(...
                Rtrain,Xtrain,Qtrain.T,Rtest,Xtest,SStot);
        otherwise
            error('bad decoder name -- jgm');
    end
    SNR(iDecoder,:) = cat(2,gather(Rsq2SNR(this_Rsqs)),NaN);
    yr_decoders(iDecoder).lagCoD =...
        getLagRsqs(this_Xhat,Rtest,Nslag,binwidth);
                
    
    % plot
    Xhat_and_t = cat(2, this_Xhat(bins,:), t');
    extra_plot_args = {yr_decoders(iDecoder).linestyle,...
        'Color',yr_decoders(iDecoder).color,'linewidth',1.5};
    for iPlot = 1:Nplots
        figure(iPlot+base_fignum);
        plot(...
            Xhat_and_t(:,kinemat_plot_inds(iPlot,1)),...
            Xhat_and_t(:,kinemat_plot_inds(iPlot,2)),...
            extra_plot_args{:});
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function finalize_reconstruction_plots(SNR, kinemat_plot_inds, kinemat,...
    yr_decoders, figure_name_prefix, figure_handles, base_fignum)
% plot true kinematics; axis labels etc.; save as tikz

%%%% TO DO:
% legend placement?


% Ns
Ndecoders = length(yr_decoders);
[Nplots,Naxes] = size(kinemat_plot_inds);
Nstates = length(kinemat);


for iPlot = 1:Nplots
    
    legend_SNR = cell(Ndecoders,1);
    figure(iPlot+base_fignum);
    
    for iAxis = 1:Naxes
        this_ind = kinemat_plot_inds(iPlot,iAxis);
        if this_ind == Nstates+1
            
            %%% assumes t is always the first variable (i.e. x axis)
            xlabel(sprintf('time (s)'))
            x_tikz_plot_str = 't';
            [legend_SNR(:)] = {'('};
            
        else
            switch kinemat(this_ind).legendname
                case 'Position'
                    units = 'mm';
                case 'Velocity'
                    units = 'mm/s';
                case 'Acceleration'
                    units = 'mm/s$^2$';
            end
            
            % collect useful elements for plot labeling, saving, etc.
            label = sprintf('%s (%s)',kinemat(this_ind).texname,units);
            tikz_plot_str = kinemat(this_ind).name;
            for iDecoder = 1:Ndecoders
                if strcmp(yr_decoders(iDecoder).name,'groundtruth')
                    legend_SNR{iDecoder} = '';
                else
                    switch iAxis
                        case 1
                            xlabel(label);
                            x_tikz_plot_str = tikz_plot_str;
                            legend_SNR{iDecoder} = sprintf('(%2.02f dB, ',...
                                SNR(iDecoder,this_ind));
                        case Naxes
                            ylabel(label);
                            y_tikz_plot_str = tikz_plot_str;
                            legend_SNR{iDecoder} = sprintf('%s%2.02f dB)',...
                                legend_SNR{iDecoder}, SNR(iDecoder,this_ind));
                        otherwise
                            error('unexpected axis number -- jgm');
                    end
                end
            end
        end
    end
    
    
    
    % now add legends and save as tikz file
    decoder_texnames = {yr_decoders.texname};
    legend_array = arrayfun...
        (@(ii)(['\tiny ', decoder_texnames{ii}, ' ', legend_SNR{ii}]),...
        1:Ndecoders,'UniformOutput',false);
    legend(legend_array,'Location','SouthEast')
    hold off
    
    [~,machine] = system('hostname');
    machine = strtrim(machine);
    switch machine
        case {'kobayashi-maru','CUPCAKE','Themistocles'}
            for iFigure = 1:length(figure_handles)
                clrNameClrPairs = arrayfun(@(ii)({yr_decoders(ii).colorname,...
                    yr_decoders(ii).color}),1:Ndecoders,...
                    'UniformOutput',false);
                matlab2tikzWrapper(sprintf(['%s_%s_vs_%s'],...
                    figure_name_prefix,y_tikz_plot_str,x_tikz_plot_str),...
                    figure_handles(iPlot),'extraColors',clrNameClrPairs);
                %'extraTikzpictureOptions',...
                %{'baseline','trim axis left','trim axis right'});
            end
        otherwise
            fprintf('skipping tikz plots on this machine -- jgm\n');
    end
    
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function model_filename = filenameForFracneurons(model_dir,...
    model_filename,fracneuron,obsv_dstrb,monkey)

% the training times are reflected in the number of hidden units/states
matching_filenames = ls([model_dir, model_filename]);
Nfiles = size(matching_filenames,1);

% find the number of hiddens in all the otherwise matching file names
numhid = zeros(Nfiles,1);
for iFile = 1:Nfiles
    this_file = matching_filenames(iFile,:);
    switch this_file(1:3)
        case 'wts'
            i0 = length(sprintf('wts_spikecounts_%s_',obsv_dstrb)) + 1;
            iF = strfind(this_file,'Hid') - 1;
        case 'LDS'
            i0 = 7;
            iF = 9;
        otherwise
            fprintf('WARNING: unexpectedly found a file name ')
            fprintf(' without LDS or wts\n');
    end
    numhid(iFile) = str2double(this_file(i0:iF));
end
switch Nfiles
    case 4 % all trainingtimes files are present
        
        % load the canonical fracneurons sweep
        load(sprintf('%sRBMish%sBMI%sRsqs_FracneuronSweep_%s_%s.mat',...
            getdir('data'),filesep,filesep,obsv_dstrb,monkey),'fracneurons')
        
        % sort Nhiddens, grab file that matches fracneuron
        [~,sorted_inds] = sort(numhid);
        iFile = sorted_inds(fracneurons==fracneuron);
        model_filename = matching_filenames(iFile,:);
        
    case 1 % probably downloaded for this purpose
        fprintf('WARNING: only found refh wts for one fracneuron.')
        fprintf('  Assuming this the correct file and proceeding...\n')
        model_filename = matching_filenames(Nfiles,:);
        
    otherwise
        fprintf('%i possible matches (out of four); don''t know',Nfiles)
        fprintf(' which fracneuron is intended so exiting\n')
        model_filename = [];
        return
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function summaryStats(decoders,kinemat,swept_param,NdataTest,...
    standardname,rivalnames,alp,swept_unit)

% Ns
Ndims = 2;
Nstates = length(kinemat);
Norder = Nstates/Ndims;
Nsessions = size(NdataTest,1);
Nswept = length(swept_param);

% indices
iStandard = strcmp({decoders(:).name},standardname);
rivals = find(ismember({decoders(:).name},rivalnames));



% table labels
row_labels = cat(2,num2str([swept_param(:)]),...
    repmat([swept_unit,' '],[length(swept_param),1]) );
row_labels = strtrim(vect(row_labels')');
col_labels = sprintf('%s ', kinemat.name);


% for each "rival"
kinnames = reshape({kinemat.legendname},Ndims,[]);
for iRival = 1:length(rivals)
    iDecoder = rivals(iRival);
    
    % SNR
    % all kinematic vars
    [MeanImprovementSNR,ISSIG] = signify(...
        Rsq2SNR(decoders(iStandard).CoD),Rsq2SNR(decoders(iDecoder).CoD),...
        NdataTest,alp,0);
    printmatJGM(squeeze(MeanImprovementSNR)',...
        sprintf('Improvement (%s over %s) in SNR (dB)',...
        decoders(iStandard).texname,decoders(iDecoder).texname),...
        row_labels,col_labels);
    printmatJGM(squeeze(ISSIG)',...
        sprintf('Significant difference (%s=:+1,%s=:-1)?',...
        decoders(iStandard).texname,decoders(iDecoder).texname),...
        row_labels,col_labels);
    
    % pooled within the plane
    [MeanImprovementSNR,ISSIG,pVal] = signify(...
        reshape(Rsq2SNR(decoders(iStandard).CoD),Nsessions*Ndims,Norder,Nswept,[]),...
        reshape(Rsq2SNR(decoders(iDecoder).CoD),Nsessions*Ndims,Norder,Nswept,[]),...
        repmat(NdataTest,[Ndims,1,1,1]),alp,0);
    printmatJGM(squeeze(MeanImprovementSNR)',...
        sprintf('Improvement (%s over %s) in SNR (dB)',...
        decoders(iStandard).texname,decoders(iDecoder).texname),...
        row_labels,sprintf('%s ',kinnames{1,:}));
    printmatJGM(squeeze(ISSIG)',...
        sprintf('Significant difference (%s=:+1,%s=:-1)?',...
        decoders(iStandard).texname,decoders(iDecoder).texname),...
        row_labels,sprintf('%s ',kinnames{1,:}));
    printmatJGM(squeeze(pVal)','Significance levels',...
        row_labels,sprintf('%s ',kinnames{1,:}));
    
    
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [MeanImprovement,ISSIG,pVal] =...
    signify(fitsStandard,fitsRival,NdataTest,alp,USEPERM)


% Ns
[Nsessions,Nstates,Nbins] = size(fitsStandard);
Nperms = 100000;

% mean differences
DiffWithRival = fitsStandard - fitsRival;
MeanImprovement = sum(NdataTest.*DiffWithRival,1)./sum(NdataTest,1);

% significance tests
if ~USEPERM
    
    % bootstrap
    Nresamples = 100000;
    
    % resample
    boot_inds = ceil(Nsessions*rand(Nsessions*Nresamples,1));
    resampled_diffs = DiffWithRival(boot_inds,:,:);
    resampled_wts   = NdataTest(boot_inds,:,:);
    
    % reshape
    resampled_diffs = permute(reshape(resampled_diffs,...
        [Nresamples,Nsessions,Nstates,Nbins]),[2,3,4,1]);
    resampled_wts   = permute(reshape(resampled_wts,...
        [Nresamples,Nsessions,1,Nbins]),[2,3,4,1]);
    
    avg_resampled_diffs = sum(resampled_wts.*resampled_diffs)./...
        sum(resampled_wts,1);
    pVal = mean(avg_resampled_diffs < 0,4);
    
    
    % based on the assumption of a normal distribution
    %     VrncImprovement = sum(NdataTest.*...
    %         (DiffWithRival - MeanImprovement).^2)./sum(NdataTest);
    %     StdErrorOfMeanImprovementSNR = sqrt(VrncImprovement/Nsessions);
    %     ISSIG = (norminv(alp/2,...
    %         abs(MeanImprovement),StdErrorOfMeanImprovementSNR) > 0).*...
    %         sign(MeanImprovement);
else
    
    % based on a permutation test
    bothFits = cat(1,fitsRival,fitsStandard);
    MdataTest = repmat(NdataTest,[2,1,1]);
    wts = permute(MdataTest./sum(MdataTest),[1,3,2]);
    inds = categorsmpl(wts(:,1),Nsessions*2*Nperms,'IndexBased');
    %%% you could decouple across bins, but it's not really necessary
    permFits = reshape(bothFits(inds,:,:),Nsessions,2,Nperms,Nstates,Nbins);
    permDstrb = permute(mean(diff(permFits,[],2)),[3,4,5,1,2]);
    pVal = mean(MeanImprovement < permDstrb);
    
end
% +1 => standard > rival, -1 => rival > standard; 0 => push
ISSIG = (pVal < alp) - (pVal > (1-alp));


% wtdFitsStandard = fitsStandard.*NdataTest./sum(NdataTest);
% wtdFitsRival = fitsRival.*NdataTest./sum(NdataTest);
%
% % cat rival first because you use diff below
% bothFits = cat(1,wtdFitsRival,wtdFitsStandard);
%
% now sample
% inds = ceil(2*Nsessions*rand(Nsessions*2*Nperms,1));
% bootFits = reshape(bothFits(inds,:,:),Nsessions,2,Nperms,Nstates,Nbins);
% bootDstrb = permute(sum(diff(bootFits,[],2)),[3,4,5,1,2]);
% pBoot = mean(MeanImprovement < bootDstrb);
% ISSIGperm = pBoot < alp;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function allDecodersBarPlots(decoders,kinemat,swept_params,NdataTest,plotversion)
% All decoders summary bar plot

% all results together
Nswept = length(swept_params);
Nstates = length(kinemat);
allRsqs = cat(4,decoders.CoD);
iStatic = strcmp({decoders(:).name},'static');


% R^2
DiffWithStatic = allRsqs(:,:,:,~iStatic) - allRsqs(:,:,:,iStatic);
MeanStaticFits = sum(NdataTest.*allRsqs(:,:,:,iStatic),1)./sum(NdataTest,1);
for iSwept = 1:Nswept
    for iState = 1:Nstates
        xaxislabels(iSwept,iState).texname =...
            sprintf('%s {\\footnotesize (%2.02f)}',...
            kinemat(iState).texname,MeanStaticFits(1,iState,iSwept));
    end
end
barPlotCore(DiffWithStatic,xaxislabels,swept_params,NdataTest,...
    'R$^2$ improvement',{decoders(~iStatic).barcolor},'CoD',plotversion,...
    'extraTikzText','[baseline,trim axis left,trim axis right]');

% SNR
DiffWithStatic = Rsq2SNR(allRsqs(:,:,:,~iStatic)) -...
    Rsq2SNR(allRsqs(:,:,:,iStatic));
MeanStaticFits = sum(NdataTest.*Rsq2SNR(allRsqs(:,:,:,iStatic)),1)./...
    sum(NdataTest,1);
for iSwept = 1:Nswept
    for iState = 1:Nstates
        xaxislabels(iSwept,iState).texname =...
            sprintf('%s {\\footnotesize (%2.02f dB)}',...
            kinemat(iState).texname,MeanStaticFits(1,iState,iSwept));
    end
end
barPlotCore(DiffWithStatic,xaxislabels,swept_params,NdataTest,...
    'SNR (dB) improvement',{decoders(~iStatic).barcolor},'SNR',plotversion,...
    'extraTikzText','[baseline,trim axis left,trim axis right]');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function barPlotCore(fitscores,xaxislabels,binwidths,NdataTest,...
    fitmetric,barcolors,fileprefix,plotversion,varargin)

% Ns
[Nsessions,Nstates,Nbinwidths,Nfilters] = size(fitscores);

% plot the averages, with error bars
MeanFits = sum(NdataTest.*fitscores,1)./sum(NdataTest,1);
VrncFits = sum(NdataTest.*(fitscores - MeanFits).^2,1)./sum(NdataTest,1);
for iBinwidth = 1:Nbinwidths
    muDiff = permute(MeanFits(1,:,iBinwidth,:),[2,4,1,3]);
    stdErr = sqrt(VrncFits(1,:,iBinwidth,:)/Nsessions);
    stdErr = permute(stdErr,[2,1,4,3]);
    
    this_xaxislabel = {xaxislabels(iBinwidth,:).texname};
    ymin = defaulter('ymin', min(0,min(vect(muDiff - stdErr))),varargin{:});
    tikzBarGraph(...
        1:Nstates,...
        muDiff,...
        cat(2,-stdErr,stdErr),...
        ymin,...
        this_xaxislabel,'kinematic variable',...
        sprintf('%s',fitmetric),'',...
        repmat(barcolors',[1,Nstates])',...
        2,0.32,'grouped',{},...
        sprintf('SLMD/%s_allvars_%03d_%s_bar',...
        fileprefix,binwidths(iBinwidth),plotversion),varargin{:});
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfBinwidthPlots(decoders,kinemat,binwidths,NdataTest,...
    SEMILOGPLOT,plotversion,rsq_lims,snr_lims)
% Within the .tex file, you may want to:
%   (1) bold the x-axis label 64 (for posters)

% Ns
Ndims = 2;                                  % two-dimensional space

% reshape
allRsqs = cat(4,decoders.CoD);
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,binwidths,NdataTest,allRsqs,'CoD',...
    rsq_lims(1),rsq_lims(2),titleStrs,filesuffices,'binwidth',...
    'bin width (ms)','R$^2$',233,SEMILOGPLOT,plotversion);

% SNR
effectOfCore(decoders,binwidths,NdataTest,Rsq2SNR(allRsqs),'SNR',...
    snr_lims(1),snr_lims(2),titleStrs,filesuffices,'binwidth',...
    'bin width (ms)','SNR (dB)',333,SEMILOGPLOT,plotversion);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfNumberOfTrainingSamples(decoders,monkeys,SEMILOGPLOT,...
    plotversion,rsq_lims,snr_lims)

% load data
[decoders,kinemat,~,trainingtimes,NdataTest] =...
    assembleData(decoders,monkeys,'Rsqs_TrainingtimeSweep_Poisson');
%assembleData(decoders,monkeys,'Rsqs_064msBins_allTrainingTimes');

summaryStats(decoders,kinemat,trainingtimes,NdataTest,'refhdynamic',...
    {'kfobs','ukf'},0.01,'s');

% reshape
Ndims = 2;
allRsqs = cat(4,decoders.CoD);
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,trainingtimes,NdataTest,allRsqs,'CoD',...
    rsq_lims(1),rsq_lims(2),titleStrs,filesuffices,'trainingtime',...
    'seconds of training data','R$^2$',433,SEMILOGPLOT,plotversion);

% SNR
effectOfCore(decoders,trainingtimes,NdataTest,Rsq2SNR(allRsqs),'SNR',...
    snr_lims(1),snr_lims(2),titleStrs,filesuffices,'trainingtime',...
    'seconds of training data','SNR (dB)',533,SEMILOGPLOT,plotversion);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfNumberOfNeurons(decoders,kinemat,binwidths,NdataTest,...
    Nneurons,mkinds,monkeys,SEMILOGPLOT,plotversion,rsq_lims,snr_lims)

% useful functions
centerByMonkey = @(XX)(cell2mat(arrayfun(@(ii)(...
    (XX(mkinds{ii},:,:,:) - mean(XX(mkinds{ii},:,:,:)))),...
    1:length(mkinds),'UniformOutput',0)'));
wtd_corrs = @(aa,bb,wt)(sum(wt.*aa.*bb)./sqrt(sum(wt.*aa.*aa).*sum(wt.*bb.*bb)));


% first see if Nneurons correlated with performance in the canonical data
allRsqs = cat(4,decoders.CoD);
[Nsessions,Nstates,~,Nfilters] = size(allRsqs);

% you have to center the data before they're comparable across animals
NneuronsCtrd = centerByMonkey(Nneurons);
reconQualityCtrd = centerByMonkey(allRsqs);

yrcorrs = reshape(wtd_corrs(reconQualityCtrd(:,:,binwidths==64,:),...
    NneuronsCtrd,NdataTest(:,:,binwidths==64)),[Nstates,Nfilters]);
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
    sprintf('SLMD/%s_allvars_%03d_%s_bar',...
    'corr_CoD_Nneurons',binwidths(3),plotversion));




reconQualityCtrd = centerByMonkey(Rsq2SNR(allRsqs));
yrcorrs = reshape(wtd_corrs(reconQualityCtrd(:,:,binwidths==64,:),...
    NneuronsCtrd,NdataTest(:,:,binwidths==64)),[Nstates,Nfilters]);
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
    sprintf('SLMD/%s_allvars_%03d_%s_bar',...
    'corr_SNR_Nneurons',binwidths(3),plotversion));



% load data
[decoders,kinemat,binwidths,trainingtimes,NdataTest,~,~,fracneurons] =...
    assembleData(decoders,monkeys,'Rsqs_FracneuronSweep_Poisson');
summaryStats(decoders,kinemat,fracneurons,NdataTest,'refhdynamic',...
    {'kfobs','ukf','wf','kfemstatic'},0.01,'units');

% reshape
Ndims = 2;
allRsqs = cat(4,decoders.CoD);
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% R^2
effectOfCore(decoders,fracneurons,NdataTest,allRsqs,'CoD',...
    rsq_lims(1),rsq_lims(2),titleStrs,filesuffices,'numneurons',...
    'fraction of all neurons','R$^2$',433,SEMILOGPLOT,plotversion);

% SNR
effectOfCore(decoders,fracneurons,NdataTest,Rsq2SNR(allRsqs),'SNR',...
    snr_lims(1),snr_lims(2),titleStrs,filesuffices,'numneurons',...
    'fraction of all neurons','SNR (dB)',533,SEMILOGPLOT,plotversion);
%%%
% At the moment you have to add this manually to the axis properties of the
% three tikz pictures:
%   xticklabel style={/pgf/number format/frac},
%%%

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfCore(decoders,xdata,NdataTest,fitscores,...
    fitmetric,fitmetricMin,fitmetricMax,titleStrs,filesuffices,sweptparam,...
    yrxlabel,yrylabel,basefignum,SEMILOGPLOT,plotversion)


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
        p = plot(xdataplot,MeanFits(iOrder,:,iDecoder));
        p.LineStyle     = decoders(iDecoder).linestyle;
        p.Color         = decoders(iDecoder).color;
        p.LineWidth     = 1.5;
        p.Marker        = decoders(iDecoder).mark;
        p.MarkerSize    = 6; % set externally? m2t uses different size sch.
        p.MarkerEdgeColor = 'none';
        p.MarkerFaceColor = 'none'; % [1.0, 1.0, 1.0]; % decoders(iDecoder).color;
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
    
    % make sure the opacity of the marks is settable
    mark_opacity_setter = '\providecommand{\theseMarkOpacities}{1}%';
    mark_code = ['\pgfplotsset{%',newline,...
        char(9),'every axis plot post/.append style={%',newline,...
        char(9),'every mark/.append style={opacity=\theseMarkOpacities}%',newline,...
        char(9),'}%',newline,'}%'];
    switch machine
        case {'kobayashi-maru','CUPCAKE','Themistocles'}
            matlab2tikzWrapper(...
                sprintf('SLMD/%s_%s_vs_%s_%s_line',...
                fitmetric,filesuffix,sweptparam,plotversion),h,...
                'extraColors',clrNameClrPairs,...
                'extraTikzpictureOptions',...
                {'baseline','trim axis left','trim axis right'},...
                'extraCode',{mark_opacity_setter, mark_code},...
                'extraAxisOptions','xticklabel style={/pgf/number format/frac}');
        otherwise
            fprintf('skipping tikz plots on this machine -- jgm\n');
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfREFHcomponents(monkeys)

plotversion = 'variantREFHs';
filename = 'Rsqs_BinwidthSweep_Poisson';

% the UKF--a baseline
decoderUKF.name = 'ukf';
decoderUKF = assembleData(decoderUKF,monkeys,filename);

% re-collect the REFH (Poisson) decoders
[decoders(1:6).name] = deal(...
    'refhstatic_hidsonly','refhdynamic_hidsonly',...
    'refhstatic_stdnrml','refhdynamic_stdnrml',...
    'refhstatic','refhdynamic');
[decoders,kinemat,binwidths,~,NdataTest] = ....
    assembleData(decoders,monkeys,filename);

% assign colors
[decoders(:).barcolor] = deal('black');
[decoders(endsWith({decoders.name},'stdnrml')).barcolor] = deal('black!50!white');
[decoders(endsWith({decoders.name},'hidsonly')).barcolor] = deal('black!20!white');
staticdecoders = {'refhstatic_hidsonly','refhstatic_stdnrml','refhstatic'};
for iDecoder = 1:length(staticdecoders)
    staticdecoder = staticdecoders{iDecoder};
    decoders = cleverAssign(decoders,staticdecoder,'barcolor',...
        [decoders(strcmp({decoders(:).name},staticdecoder)).barcolor,...
        ',postaction={pattern=north east lines,pattern color=white}']);
end
allRsqs = cat(4,decoders.CoD);




% Ns
[Nsessions, Nstates, Nbinwidths, Ndecoders] = size(allRsqs);
Ndims = 2; % fact

% reshape to pool across dimensions, within variables (pos,vel,acc)
NdataTest = repmat(NdataTest,[Ndims,1,1]);
doubleup_func = @(TT)(reshape(TT,Nsessions*Ndims,Nstates/Ndims,...
    Nbinwidths,[]));
rsqs_all = doubleup_func(allRsqs);
rsqs_ukf = doubleup_func(decoderUKF.CoD);



% R^2
diff_with_ukf = rsqs_all - rsqs_ukf;
MeanStaticFits = sum(NdataTest.*rsqs_ukf,1)./sum(NdataTest,1);
for iSwept = 1:Nbinwidths
    for iState = Ndims:Ndims:Nstates
        xaxislabels(iSwept,iState/Ndims).texname =...
            sprintf('%s {\\footnotesize (%2.02f)}',...
            kinemat(iState).legendname,...
            MeanStaticFits(1,iState/Ndims,iSwept));
    end
end
barPlotCore(diff_with_ukf,xaxislabels,binwidths,NdataTest,...
    'R$^2_\text{rEFH}$ - R$^2_\text{UKF}$',...
    {decoders.barcolor},'CoD_rEFH',plotversion);


% SNR
bruteforcetext = getbruteforcetext;
diff_with_ukf = Rsq2SNR(rsqs_all) - Rsq2SNR(rsqs_ukf);
MeanStaticFits = sum(NdataTest.*Rsq2SNR(rsqs_ukf),1)./sum(NdataTest,1);
for iSwept = 1:Nbinwidths
    for iState = Ndims:Ndims:Nstates
        xaxislabels(iSwept,iState/Ndims).texname =...
            sprintf('%s (%2.02f dB)',...
            kinemat(iState).legendname,...
            MeanStaticFits(1,iState/Ndims,iSwept));
    end
end
barPlotCore(diff_with_ukf,xaxislabels,binwidths,NdataTest,...
    'SNR$_\text{rEFH}$ - SNR$_\text{UKF}$ (dB)',...
    {decoders.barcolor},'SNR_rEFH',plotversion,...
    'extraAxisText',bruteforcetext,'ymin',-0.5);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tappedDecoding(monkeys)

% untapped and tapped versions
[decoders(1:4).name] = deal('ukf','refhdynamic',...
    'ukf321','refhdynamic320_hidsonly');

% load the saved results
[decoders,kinemat,binwidths,trainingtimes,NdataTest,Nneurons,mkinds] =...
    assembleData(decoders,monkeys,'Rsqs_BinwidthSweep_Poisson');
fprintf('"canonical" training time is %d seconds\n',trainingtimes);

% for the rEFH don't use tapped decoding beyond 64 ms
decoders(strcmp({decoders(:).name},'refhdynamic320_hidsonly')).CoD(...
    :,:,binwidths > 64) =...
    decoders(strcmp({decoders(:).name},'refhdynamic')).CoD(...
    :,:,binwidths > 64);

% plot the effect of bin widths
effectOfBinwidthPlots(decoders,kinemat,binwidths,NdataTest,0,...
    'tappeddecoders',[0,0.8],[0,6]);

% significant differences??
summaryStats(decoders,kinemat,binwidths,NdataTest,...
    'refhdynamic320_hidsonly',{'ukf321'},0.01,'ms');

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
decoders = cleverAssign(decoders,'wf',          'texname','WF');
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

decoders = cleverAssign(decoders,'refhdynamic320_hidsonly','texname',...
    'rEFH (+KF, tapped)');
decoders = cleverAssign(decoders,'ukf321','texname','UKF (tapped)');

% ground truth
decoders = cleverAssign(decoders,'groundtruth','texname','actual');


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignColorNames(decoders)
% color names (order of these statements doesn't matter)
%
% NB that this is actually more important than assignColors, since the
% color name written into the tex file by this macro will ultimately be
% defined externally in a master tex file, overriding the color definitions
% written on the basis of setColor.m.

% standard decoders
decoders = cleverAssign(decoders,'static',     'colorname','xtraclr');
decoders = cleverAssign(decoders,'kfobs',      'colorname','OBSclr');
decoders = cleverAssign(decoders,'kfemstatic', 'colorname','EMclr');
decoders = cleverAssign(decoders,'kfemdynamic','colorname','EMclr');
decoders = cleverAssign(decoders,'wf',         'colorname','mustardclr');
decoders = cleverAssign(decoders,'ukf',        'colorname','pinkish');
decoders = cleverAssign(decoders,'refhstatic', 'colorname','EFHclr');
decoders = cleverAssign(decoders,'refhdynamic','colorname','EFHclr');

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml','colorname','EFHclr');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml','colorname','EFHclr');
decoders = cleverAssign(decoders,'refhstatic_hidsonly','colorname','gray');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','colorname','gray');

decoders = cleverAssign(decoders,'refhdynamic320_hidsonly','colorname','EFHclr');
decoders = cleverAssign(decoders,'ukf321','colorname','pinkish');

% ground truth
decoders = cleverAssign(decoders,'groundtruth','colorname','black');


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
decoders = cleverAssign(decoders,'wf',          'color',mustardcolor);
decoders = cleverAssign(decoders,'ukf',         'color',pinkish);
decoders = cleverAssign(decoders,'refhstatic',  'color',EFHcolor);
decoders = cleverAssign(decoders,'refhdynamic', 'color',EFHcolor);


% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml',  'color',EFHcolor);
decoders = cleverAssign(decoders,'refhdynamic_stdnrml', 'color',EFHcolor);
decoders = cleverAssign(decoders,'refhstatic_hidsonly', 'color',[0.2,0.2,0.2]);
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','color',[0.2,0.2,0.2]);

decoders = cleverAssign(decoders,'refhdynamic320_hidsonly','color',EFHcolor);
decoders = cleverAssign(decoders,'ukf321','color',pinkish);

% ground truth
decoders = cleverAssign(decoders,'groundtruth','color',[0,0,0]);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignMarks(decoders)
% Marks (basically for a black-and-white verion of the paper)

%%% MARKS
% standard decoders
decoders = cleverAssign(decoders,'static',     'mark','x');
decoders = cleverAssign(decoders,'kfobs',      'mark','o');
decoders = cleverAssign(decoders,'kfemstatic', 'mark','s');
decoders = cleverAssign(decoders,'kfemdynamic','mark','s');
decoders = cleverAssign(decoders,'wf',         'mark','*');
decoders = cleverAssign(decoders,'ukf',        'mark','v');
decoders = cleverAssign(decoders,'refhstatic', 'mark','^');
decoders = cleverAssign(decoders,'refhdynamic','mark','^');
    
% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml',  'mark','^');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml', 'mark','^');
decoders = cleverAssign(decoders,'refhstatic_hidsonly', 'mark','^');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','mark','^');

decoders = cleverAssign(decoders,'refhdynamic320_hidsonly', 'mark','^');
decoders = cleverAssign(decoders,'ukf321',                  'mark','v');

% ground truth
decoders = cleverAssign(decoders,'groundtruth','mark','none');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignPlotLineTypes(decoders)
% the plot line types (order of these statements doesn't matter)

% standard decoders
decoders = cleverAssign(decoders,'static',      'linestyle','-');
decoders = cleverAssign(decoders,'kfobs',       'linestyle','-');
decoders = cleverAssign(decoders,'kfemstatic',  'linestyle','--');
decoders = cleverAssign(decoders,'kfemdynamic', 'linestyle','-');
decoders = cleverAssign(decoders,'wf',          'linestyle','-');
decoders = cleverAssign(decoders,'ukf',         'linestyle','-');
decoders = cleverAssign(decoders,'refhstatic',  'linestyle','--');
decoders = cleverAssign(decoders,'refhdynamic', 'linestyle','-');

% non-standard decoders
decoders = cleverAssign(decoders,'refhstatic_stdnrml','linestyle','--');
decoders = cleverAssign(decoders,'refhdynamic_stdnrml','linestyle','-');
decoders = cleverAssign(decoders,'refhstatic_hidsonly','linestyle','--');
decoders = cleverAssign(decoders,'refhdynamic_hidsonly','linestyle','-');
%%% you would perhaps like to alternate colors....

decoders = cleverAssign(decoders,'refhdynamic320_hidsonly','linestyle',':');
decoders = cleverAssign(decoders,'ukf321','linestyle',':');

% ground truth
decoders = cleverAssign(decoders,'groundtruth','linestyle',':');


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function decoders = assignBarColorStyles(decoders)
% the bar color styles (order of these statements doesn't matter)
%
% The "static" decoders are hatched.

[decoders(:).barcolor] = deal(decoders(:).colorname);

%%% you *could* do this based on the 'linestyle'---assuming
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



%%%% tapped decoders


%%%%

%%%% ground truth

%%%%



end
%-------------------------------------------------------------------------%

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







%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function allDecodersBarPlots_one_monkey(monkey)

% this sets the names but also the "canonical" ordering!
[decoders(1:8).name] = deal('static','kfobs','kfemstatic',...
    'kfemdynamic','wf','ukf','refhstatic','refhdynamic');

% load the saved results
[decoders,kinemat,binwidths,trainingtimes,NdataTest,Nneurons,mkinds] =...
    assembleData(decoders,{monkey},'Rsqs_BinwidthSweep_Poisson');
fprintf('"canonical" training time is %d seconds\n',trainingtimes);

% bar summary plots (\FigAllDecodersBarPlots)
allDecodersBarPlots(decoders,kinemat,binwidths,NdataTest,...
    sprintf('%s_%s','alldecoders',monkey));


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function compare_LGDS_with_Fixed_20_Hids
% 20 hidden units vs. Nneurons/3

clear; clc;

% useful
yr_decoder = 'kfemstatic';

% 20 hidden units
I20 = load(sprintf('%sRBMish/BMI/Rsqs_BinwidthSweep_Poisson_Indy_20hids.mat',...
    getdir('data')));
Rsqs_20 = I20.Rsqs(:,:,:,strcmp(I20.decodernames,yr_decoder));


% Nneurons/3
Ivariable = load('C:\#DATA\RBMish/BMI/Rsqs_BinwidthSweep_Poisson_Indy.mat');
Rsqs_variable = Ivariable.Rsqs(:,:,Ivariable.binwidths==I20.binwidths,...
    strcmp(Ivariable.decodernames,yr_decoder));

% plot
title_strs = {'position', 'position',...
    'velocity', 'velocity', 'acceleration', 'acceleration'};
for iState = 1:2:size(Rsqs_variable,2)
    
    hh = figure(iState); clf; hold on
    scatter(Rsqs_variable(:,iState),Rsqs_20(:,iState),...
        150,'^','MarkerEdgeColor','k','MarkerFaceColor','none');
    scatter(Rsqs_variable(:,iState+1),Rsqs_20(:,iState+1),...
        'o','MarkerEdgeColor','k','MarkerFaceColor','none');
    xlabel('R$^2$: Nhid=Nneurons/3');
    ylabel('R$^2$: Nhid=20');
    plot([0,1],[0,1],'k')
    legend('lateral','sagittal','Location','NorthWest');
    title(title_strs{iState})
    axis equal tight
    
    matlab2tikzWrapper(sprintf('SLMD/%s_Nhid_var_vs_Nhid_20_%s',...
        yr_decoder,title_strs{iState}),hh);
    
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function compare_EM4LDS_jgm_vs_zg
% EM for LGDS: JGM vs. ZG

hh = open('C:\#DATA\RBMish\BMI\xentropy_jgm_zg.fig');
xlabel('iterations of EM')
ylabel('cross entropy')
title('EM for LGDS: JGM vs. ZG');

matlab2tikzWrapper(sprintf('SLMD/kfemstatic_xentropy_Indy_Session_15_jgm_zg'),hh);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectOfNumberOfNeurons_scatter(decoders,kinemat,binwidths,...
    Nneurons,mkinds,monkeys,decodernames)

%
M1plusS1inds = 1:7;  % fact, but would be nicer to extract from raw...
allSNRs = Rsq2SNR(cat(4,decoders.CoD));
Nmsperbin = 64;

% Ns
Nstates = length(kinemat);
Nmonkeys = length(mkinds);

% cbrewer Dark2
clrs = [
    27,158,119;...      % greenish
    217,95,2;...        % orangish
    117,112,179;...     % purplish
    231,41,138;...      % magenta
    ]/255;


for iDecoder = 1:length(decodernames)
    yr_decoder = decodernames{iDecoder};
    for iState = 1:2:Nstates
        h = figure(1005+iState); clf; hold on;
        for iMk = 1:Nmonkeys
            if strcmp(monkeys(iMk),'Indy')
                yr_inds = setdiff(mkinds{iMk},M1plusS1inds);
            else
                yr_inds = mkinds{iMk};
            end
            scatter(...
                Nneurons(yr_inds),...
                allSNRs(yr_inds,iState,binwidths==Nmsperbin,...
                strcmp({decoders.name},yr_decoder)),...
                150,'^','MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
            scatter(...
                Nneurons(yr_inds),...
                allSNRs(yr_inds,iState+1,binwidths==Nmsperbin,...
                strcmp({decoders.name},yr_decoder)),...
                'o','MarkerEdgeColor',clrs(iMk,:),'MarkerFaceColor','none');
            
            %labels = arrayfun(@(ii)([monkeylabels{ii}, ', ',varsymbols{ii}]),...
            %    1:numel(monkeylabels),'UniformOutput',0);
            xlabel('number of neurons');
            ylabel(sprintf('SNR (dB), %s',...
                decoders(strcmp({decoders.name},yr_decoder)).texname))
            title(sprintf('%s (%ims bins)',kinemat(iState).legendname,Nmsperbin))
            axis square
        end
        yr_inds = M1plusS1inds;
        scatter(...
            Nneurons(yr_inds),...
            allSNRs(yr_inds,iState,binwidths==Nmsperbin,...
            strcmp({decoders.name},yr_decoder)),...
            150,'^','MarkerEdgeColor',clrs(Nmonkeys+1,:),'MarkerFaceColor','none');
        scatter(...
            Nneurons(yr_inds),...
            allSNRs(yr_inds,iState+1,binwidths==Nmsperbin,...
            strcmp({decoders.name},yr_decoder)),...
            'o','MarkerEdgeColor',clrs(Nmonkeys+1,:),'MarkerFaceColor','none');
        
        hold off;
        
        matlab2tikzWrapper(sprintf('SLMD/SNR_%s_%0.3i_%s_vs_Nneurons_scatter',...
            lower(kinemat(iState).legendname(1:3)),Nmsperbin,yr_decoder),h);
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function effectsOfSparsityTargetandNumberOfHiddenUnits(monkeys)

% load data
[decoders(1:2).name] = deal('refhstatic','refhdynamic');
[decoders,kinemat,~,~,NdataTest,~,~,~,hidsensoryratios,phidtargets] =...
    assembleData(decoders,monkeys,...
    'Rsqs_PhidtargetSweep_HidsensoryratioSweep_Poisson');

% reshape
allRsqs = cat(4,decoders.CoD);
Ndims = 2;
[Nsessions,Nstates,Nsparsitytargets,Ndecoders,Nhidsensoryratios] =...
    size(allRsqs);
allRsqs = reshape(allRsqs,[Nsessions,Ndims,Nstates/Ndims,...
    Nsparsitytargets,Ndecoders,Nhidsensoryratios]);
titleStrs = reshape({kinemat.legendname},Ndims,[]);
filesuffices = reshape({kinemat.name},Ndims,[]);

% useful function
reconQuality_func = @(RR)(reshape(mean(sum(NdataTest.*RR,1)./sum(NdataTest,1),2),...
    [1,Nstates/Ndims,Nsparsitytargets,Ndecoders,Nhidsensoryratios]));

% plot
meanCoD = reconQuality_func(allRsqs);
meanSNR = reconQuality_func(Rsq2SNR(allRsqs));

for iDecoder = 1:Ndecoders
    for iDeriv = 1:(Nstates/Ndims)
        hf = figure;
        imagesc(squeeze(meanCoD(:,iDeriv,:,iDecoder,:)));
        hc = colorbar;
        ylabel(hc, 'R$^2$')
        title(titleStrs{1,iDeriv});
        set(gca, 'XTick', 1:Nhidsensoryratios, 'XTickLabel', hidsensoryratios)
        xlabel('$N_\text{hid}/N_\text{sensory}$')
        set(gca, 'YTick', 1:Nsparsitytargets,'YTickLabel', phidtargets)
        ylabel('sparsity target')
        
        matlab2tikzWrapper(sprintf(...
            'SLMD/%s_%s_vs_phidtargets_and_hidsensoryratios_%s_heatmap',...
            'CoD',filesuffices{1,iDeriv}(1:end-1),decoders(iDecoder).name),...
            hf,'relativeDataPath','../../tikzpics/SLMD');
        
        
        hf = figure;
        imagesc(squeeze(meanSNR(:,iDeriv,:,iDecoder,:)));
        hc = colorbar;
        ylabel(hc, 'SNR (dB)')
        title(titleStrs{1,iDeriv});
        set(gca, 'XTick', 1:Nhidsensoryratios, 'XTickLabel', hidsensoryratios)
        xlabel('$N_\text{hid}/N_\text{sensory}$')
        set(gca, 'YTick', 1:Nsparsitytargets,'YTickLabel', phidtargets)
        ylabel('sparsity target')
        
        matlab2tikzWrapper(sprintf(...
            'SLMD/%s_%s_vs_phidtargets_and_hidsensoryratios_%s_heatmap',...
            'SNR',filesuffices{1,iDeriv}(1:end-1),decoders(iDecoder).name),...
            hf,'relativeDataPath','../../tikzpics/SLMD');
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function Rsqs = getLagRsqs(X,R,Nslag,Nmsperbin)

% Ns
Nmspers = 1000;
Nbinslag = floor(Nslag*Nmspers/Nmsperbin);
[Nbins,Nstates] = size(X);

% malloc
Rsqs = zeros(Nbinslag,Nstates,'like',R);

% cycle through lags
for tau = 0:(Nbinslag-1)
    Mbins = Nbins - tau;
    [~,Rsqs(tau+1,:)] = linregress(...
        [ R(1:Mbins,:) ,ones(Mbins,1)], X(tau+(1:Mbins),:));
end


% figure()
% plot( (0:-1:(1-Nbinslag))*Nmsperbin/Nmspers, Rsqs )
% ylabel('Rsq')
% xlabel('time (s)')
% ax = axis;
% axis([ax(1), ax(2), 0, 1.0]);
% title(titlestr)


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotReconstructionLaggedSpikeCoDs(decoders,kinemat,Nmsperbin)

Ndims = 2;
Nmspers = 1000;
[Nbinslag,Nstates] = size(decoders(1).lagCoD);
Nderiv = Nstates/Ndims;
t = (0:-1:(1-Nbinslag))*Nmsperbin/Nmspers;


% set up
titlemat = reshape({kinemat.legendname},Ndims,Nderiv);
figoffset = 31341;
for iDeriv=1:Nderiv, figure(figoffset+iDeriv); clf; hold on; end

% plot each decoder
for iDecoder = 1:length(decoders)
    lagCoDs = mean(reshape(decoders(iDecoder).lagCoD,...
        [Nbinslag,Ndims,Nderiv]),2);
    for iDeriv = 1:Nderiv
        figure(figoffset+iDeriv)
        plot(t, lagCoDs(:,:,iDeriv),decoders(iDecoder).linestyle,...
            'Color',decoders(iDecoder).color,'linewidth',1.5);
    end
    
end


clrNameClrPairs = arrayfun(@(ii)({decoders(ii).colorname,...
    decoders(ii).color}),1:length(decoders),...
    'UniformOutput',false);
for iDeriv=1:Nderiv
    hh = figure(figoffset+iDeriv); 
    legend({decoders.texname},'Location','NorthWest')
    ax = axis;
    axis([ax(1), ax(2), 0, 0.8])
    title(titlemat{1,iDeriv})
    ylabel('R$^2$');
    xlabel('time (s)')
    hold off;
    
    matlab2tikzWrapper(...
        sprintf('SLMD/CoD_%s_%03i_reconlaggedspike_alldecoders',...
        lower(titlemat{1,iDeriv}(1:3)),Nmsperbin),hh,...
        'extraColors',clrNameClrPairs);
end

end
%-------------------------------------------------------------------------%


















































