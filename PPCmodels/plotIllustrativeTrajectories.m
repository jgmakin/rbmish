function plotIllustrativeTrajectories(datatype,algorithm,Ndims,VARIANTS)
% plotIllustrativeTrajectories  Just like it says
%
% USAGES:
%   plotIllustrativeTrajectories('LTI-PPC','rEFH',1,{})
%   plotIllustrativeTrajectories('LTI-PPC','rEFH',1,{'withEC'})
% 

%-------------------------------------------------------------------------%
% Revised: 03/04/16
%   -functionized, rationalized
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%

[wts,params] = getModelwtsparams(algorithm,datatype,Ndims,VARIANTS);
plotRandomTrajectoryPiece(1000,300:600,wts,params);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [wts,params] = getModelwtsparams(algorithm,datatype,Ndims,VARIANTS)

variantstr = arrayfun(@(i)(['_',VARIANTS{i}]),1:length(VARIANTS),...
    'UniformOutput',false);
yrfilter = sprintf(['%s_%iD_%s',variantstr{:}],algorithm,Ndims,datatype);
load([getdir('data'),'RBMish/EFHs/wts_',yrfilter,'_ManyXprmts.mat'])
wts = Allwts{1};

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotRandomTrajectoryPiece(T,inds,wts,params)

% init
dt = params.dynamics.dt;
tvec = dt*inds;
setColors;
yticklabelArray(:,1) = {'$-2\pi/6$','$-\pi/6$','$0$','$\pi/6$','$2\pi/6$'};
yticklabelArray(:,2) = {'$-1.25$','$-0.625$','$0$','$0.625$','$1.25$'};
RGBcolor(:,1) = PROPcolor;
RGBcolor(:,2) = ECcolor;
colorNames{1} = 'pclr';
colorNames{2} = 'efcpclr';
cbrewerschemes{1} = 'Oranges';
cbrewerschemes{2} = 'Purples';
N = params.N;
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

% generate data
[X,Q] = params.getLatents(1*T,dataclass);
R = params.getData(X,Q);

% decode inputs
[Shat0,ttlSpks0] = decodeDataPPC(R,X,Q,params);
Info0 = GTPNposteriorInfo(ttlSpks0,params);
pSENSORY = cumulantNeutralize(Shat0,Info0,params);

% filter with rEFH
[~,~,Shat,Info] = testEFHPPC(R,X,Q,wts,params);
pREFH = cumulantNeutralize(Shat,Info,params);

for iMod = 1:length(params.mods)
    
    % plot in two panels
    fighandle = figure(1021+iMod); clf
    p = panel();
    p.pack(2,1);
    
    % min and max for this modality
    smin = params.smin(:,iMod);
    smax = params.smax(:,iMod);
    
    % panel 1: line plot
    p(1,1).select();
    hold on;
    plot(tvec,squeeze(pSENSORY.Xpct(inds,:,iMod))',...
        'Color',RGBcolor(:,iMod),'Linewidth',1.5)
    plot(tvec,squeeze(pREFH.Xpct(inds,:,iMod))',...
        'Color',EFHcolor,'Linewidth',1.5)
    plot(tvec,squeeze(X(inds,iMod*2-1))',...
        'Color',[0.4 0.4 0.4],'Linewidth',1.5)
    % plot(tvec,squeeze(pKFtrue.Xpct(1,1,:))','color',OPTcolor)
    xticks = linspace(tvec(1),tvec(end),4);
    set(gca,'XTick',xticks)
    ylabel('\prop (rad)','Interpreter','none')
    yticks = linspace(smin,smax,5);
    set(gca,'YTick',yticks)
    set(gca,'YTickLabel',yticklabelArray(:,iMod));
    axis([tvec(1),tvec(end),smin,smax])
    hold off;
    
    % panel 2: heatmap
    p(2,1).select();
    colormap(cbrewer('seq',cbrewerschemes{iMod},max(R(:))));
    %%%% make sure min = white!
    image(R(inds,(1:N)+N*(iMod-1))');
    axis xy; axis tight;
    xlabel('time (s)')
    set(gca,'XTickLabel','')
    ylabel('neuron no.')
    set(gca,'YTick',[1,5,10,15]); %%% hard coded!!
    colorbar('East');
    
    
    % wrap into tikz
    figname = ['illustrativeTraj_',params.datatype,'_',params.mods{iMod},'_',date];    
    matlab2tikzWrapper(figname,fighandle,'extraColors',...
        { {colorNames{iMod},RGBcolor(:,iMod)}, {'rbmclr',EFHcolor} });
    %%% you still have to:
    % (1) change the directory to the .png file
    %       --but this can be fixed by using:
    %           matlab2tikz('dataPath',CHAR, ...)
    % (2) for *both* subplots, change something like:
    %   at={(0.590551in,0.590551in)},
    % to something like
    %   at={(0.590551in,\thisTikzPicScale*0.590551in)},
    % (3) colorbar.....
end


end
%-------------------------------------------------------------------------%