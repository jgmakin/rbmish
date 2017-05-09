function illustrateUnisensoryResponses(r,params)
% This function draws some of the illustrative subfigures for yr paper.
% It produces figures very similar to drawnoisyhills.m's, but it puts
% tuning curves on the plots.
%
% USAGE:
%   load([getdir('data'),'RBMish/EFHs/wts_2Dinteg_160226.mat'])
%   R = generateData(1000,params.getLatents,params.getData,'double');
%   illustrateUnisensoryResponses(R(ceil(rand*size(R,1)),:),params)

%-------------------------------------------------------------------------%
% Revised: 08/24/16
%   -changed to accommodate grid2world's newly transposed input and output.
%   -replaced reference to decode.m with GTPNsuffstats.m
% Revised: 02/29/16
%   -massively rewrote: now writes a tikz file directly, rather than using
%   matlab2tikz
% Revised: 02/26/16
%   -improved; in particular, now uses colors from getColor.m rather than
%   approximations
% Cribbed: 12/03/14 (happy b'day, DN)
%   -from drawnoisyhills.m
%-------------------------------------------------------------------------%

% init
Nmods = length(params.mods);
smin = params.smin;
smax = params.smax;
N = params.N;
tuningCov = computetuningcovs(params);

% transform...
cntrsOfMass = GTPNsuffstats(r,params);
T = displayshape(r,params);

% loop through modalities
for iMod = 1:Nmods
    
    % stimuli, example responses (as a grid), and example center of mass
    S = grid2world([1:N;1:N]',[smin(:,iMod),smax(:,iMod)],params);
    Rmod = T{iMod};
    
    % color map for this modality
    clr = getColor(params.mods{iMod});
    clrmap = [...
        linspace(1,clr(1),max(Rmod(:)))',...
        linspace(1,clr(2),max(Rmod(:)))',...
        linspace(1,clr(3),max(Rmod(:)))'];
    
    
    % different outlines for different modalities
    switch params.mods{iMod}
        case 'Hand-Position'
            scale = 0.1;    
            
            extraAxisOptionsStr = [...
                'xlabel={$\visl_1$ (transverse hand pos., cm)},',sprintf('\n\t'),...
                'ylabel={$\visl_2$ (sagittal hand pos., cm)},',sprintf('\n\t')];
            outlineStr = tikzReachableWorkspace(params.roboparams);
            tikzfilename = ['exVislResp-',date];
            
        case 'Joint-Angle'
            scale = 2.5;
            
            % hard-coded :-(
            extraAxisOptionsStr = [...
                'xtick={-1.5707963267949,-0.785398163397448,0,0.785398163397448},',sprintf('\n\t'),...
                'xticklabels={{-$\pi$/2},{-$\pi$/4},{0},{$\pi$/4}},',sprintf('\n\t'),...
                'xlabel={$\prop_1$ (shoulder angle, rad)},',sprintf('\n\t'),...
                'ytick={0.785398163397448,1.5707963267949,2.35619449019234},',sprintf('\n\t'),...
                'yticklabels={{$\pi$/4},{$\pi$/2},{3$\pi$/4}},',sprintf('\n\t'),...
                'ylabel={$\prop_2$ (elbow angle, rad)},',sprintf('\n\t')];
            outlineStr = ['\draw (',...
                num2str(params.roboparams.thmin(1)),',',...
                num2str(params.roboparams.thmin(2)),') rectangle (',...
                num2str(params.roboparams.thmax(1)),',',...
                num2str(params.roboparams.thmax(2)),');',sprintf('\n')];
            tikzfilename = ['exPropResp-',date];
            
        otherwise
            error('you never coded up this case -- jgm\n');
    end
    
    % (1) the heat map of responses
    outtxt = tikzHeatMap(S(:,1),S(:,2),Rmod,clrmap,scale,extraAxisOptionsStr);
    
    % (2) the center of mass
    outtxt = outtxt(1:(end-length('\end{axis}')-1)); 
    outtxt = [outtxt,sprintf('\n'),...
        '\addplot [color=black,mark size=4.0pt,only marks,mark=x,'...
        'mark options={solid},forget plot]',sprintf('\n'),...
        'table[row sep=crcr]{%',sprintf('\n'),...
        num2str(cntrsOfMass(1,1,iMod)),sprintf('\t'),...
        num2str(cntrsOfMass(1,2,iMod)),'\\',sprintf('\n'),...
        '};',sprintf('\n'),'\end{axis}',sprintf('\n')];
    
    % (3) the workspace
    outtxt = [outtxt,outlineStr];
    
    % (4) the tuning iso-contours (ellipses)
    [inds1,inds2] = find(Rmod>0);
    for iResp = 1:length(inds1)
        tuningStd = chol(tuningCov{iMod})';
        Cntr = [S(inds1(iResp),1) S(inds2(iResp),2)];
        outtxt = [outtxt,tikzEllipse(Cntr,tuningStd,'lightgray',0.5)];
    end
    
    % write out
    scaleStr = ['[x=',num2str(scale),'cm,y=',num2str(scale),'cm]'];
    tikzWrite(outtxt,tikzfilename,scaleStr);
   
  
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%%%% retired
function plotHeatMap(x,y,Z,clr,k)

% imagesc(x,y,Z);
% set(gca,'YDIR','normal')
x = linspace(x(1),x(2),size(Z,1));
xspacing = x(end)-x(end-1);
x = [x, x(end) + xspacing] - xspacing/2;

y = linspace(y(1),y(2),size(Z,1));
yspacing = y(end)-y(end-1);
y = [y, y(end) + yspacing] - yspacing/2;

tiny = 0.0001; %%% this accounts for a bug in matlab2tikz, you think
pcolor(x,y,Z([1:end,end],[1:end,end])+tiny);
axis xy
shading flat
%%% shading faceted

clrmap = [...
    linspace(1,clr(1),max(Z(:)))',...
    linspace(1,clr(2),max(Z(:)))',...
    linspace(1,clr(3),max(Z(:)))'];
colormap(k,clrmap)
colorbar;
axis tight


%%%%% saveas(gcf,'test','epsc');
% some colorbars won't export properly to .eps if you save from the file
% menu; this saveas command somehow overcomes that problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%-------------------------------------------------------------------------%