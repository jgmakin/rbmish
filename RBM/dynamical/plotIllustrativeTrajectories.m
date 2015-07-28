% plotIllustrativeTrajectories

load('dynamical\finalwts\wts1DrEFHManyXprmts.mat')
wts = Allwts{1};
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;
dt = params.dynamics.dt;
T = params.dynamics.T;

inds = 300:600;
tvec = dt*inds;
setColors;


LDSdata = getLDSdata(params);

%%
fighandle = figure(1021); clf
iCase = ceil(params.Ncases*rand);
p = panel();
p.pack(2,1);


p(1,1).select();
hold on;
plot(tvec,squeeze(LDSdata.Z(iCase,1,inds))','Color',[0.4 0.4 0.4],'Linewidth',1.5)
plot(tvec,squeeze(LDSdata.Y(iCase,1,inds))','Color',PROPcolor,'Linewidth',1.5)
% plot(tvec,squeeze(pKFtrue.Xpct(iCase,1,:))','color',OPTcolor)
xticks = linspace(tvec(1),tvec(end),4);
set(gca,'XTick',xticks)
ylabel('\prop (rad)','Interpreter','none')
yticks = params.thmin:pi/6:params.thmax;
set(gca,'YTick',yticks)
set(gca,'YTickLabel',{'$-2\pi/6$','$-\pi/6$','$0$','$\pi/6$','$2\pi/6$'});
axis([tvec(1),tvec(end),params.thmin,params.thmax])
hold off;


p(2,1).select();
colormap(cbrewer('seq','Oranges',max(LDSdata.R(:))))
image(squeeze(LDSdata.R(iCase,:,inds)));
axis xy; axis tight;
xlabel('time (s)')
set(gca,'XTickLabel','')
ylabel('neuron no.')
set(gca,'YTick',[1,5,10,15]); %%% hard coded!!
colorbar('East');

%%
matlab2tikzWrapper(['illustrativeTraj',params.MODEL,date],fighandle);
%%% you still have to:
% (1) change the directory to the .png file
% (2) tighten up the axis labels:
%   ylabel absolute, ylabel style={yshift=-0.5em},







%% CONTROLLED NETWORK
clear; clc; 

load('dynamical\finalwts\wts1DrEFHwithECManyXprmts.mat');
wts = Allwts{1};
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;
dt = params.dynamics.dt;
T = params.dynamics.T;

inds = 600:900;
tvec = dt*inds;
setColors;


LDSdata = getLDSdata(params);




%% 
iCase = ceil(params.Ncases*rand);
yticklabelArray(:,1) = {'$-2\pi/6$','$-\pi/6$','$0$','$\pi/6$','$2\pi/6$'};
yticklabelArray(:,2) = {'$-1.25$','$-0.625$','$0$','$0.625$','$1.25$'};
RGBcolor(:,1) = PROPcolor;
RGBcolor(:,2) = ECcolor;
colorNames{1} = 'pclr';
colorNames{2} = 'efcpclr';
cbrewerschemes{1} = 'Oranges';
cbrewerschemes{2} = 'Purples';
N = params.N;

for iMod = 1:params.Nmods
   
    fighandle = figure(1021+iMod); clf
    pCTRL = panel();
    pCTRL.pack(2,1);
    
    smin = params.smin(:,iMod);
    smax = params.smax(:,iMod);


    % trajectories
    pCTRL(1,1).select();
    hold on;
    plot(tvec,squeeze(LDSdata.Z(iCase,iMod*2-1,inds))',...
        'Color',[0.4 0.4 0.4],'Linewidth',1.5)
    plot(tvec,squeeze(LDSdata.Y(iCase,iMod,inds))',...
        'Color',RGBcolor(:,iMod),'Linewidth',1.5)
    % plot(tvec,squeeze(pKFtrue.Xpct(iCase,1,:))','color',OPTcolor)
    xticks = linspace(tvec(1),tvec(end),4);
    set(gca,'XTick',xticks)
    ylabel('\prop (rad)','Interpreter','none')
    yticks = linspace(smin,smax,5);
    set(gca,'YTick',yticks)
    set(gca,'YTickLabel',yticklabelArray(:,iMod));
    axis([tvec(1),tvec(end),smin,smax])
    hold off;
    
    
    % neurons (heat map)
    pCTRL(2,1).select();
    colormap(cbrewer('seq',cbrewerschemes{iMod},max(LDSdata.R(:))))
    %%%% make sure min = white!
    image(squeeze(LDSdata.R(iCase,(1:N)+N*(iMod-1),inds)));
    axis xy; axis tight;
    xlabel('time (s)')
    set(gca,'XTickLabel','')
    ylabel('neuron no.')
    set(gca,'YTick',[1,5,10,15]); %%% hard coded!!
    colorbar('East');
    
    figname = ['illustrativeTraj',params.MODEL,params.mods{iMod},date];
    %%%matlab2tikzWrapper(figname,fighandle,{colorNames{iMod},RGBcolor(:,iMod)});
    %%% you still have to:
    % (1) change the directory to the .png file
    % (2) tighten up the axis labels:
    %   ylabel absolute, ylabel style={yshift=-0.5em},
    % (3) change colors to colors defined in rbmish.sty

    
end


