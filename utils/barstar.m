function h = barstar(data,xlbl,ylbl,xTickStrs,yTickStrs,...
    barClrs,titlestr,fontsize,ymax)
% A better bar plotter.
%   One often wants to make a bar plot in which
%       (1) the different bars are different colors, and
%       (2) the different bars have different labels, very often using Tex.
%   This function allows all of that to be done.
%
% Note that this involves working around an annoying issue in Matlab (or 
% really the renderer): changing the bar colors normally uses RGB colors
% for the facevertexcdata, but RGB is incompatible with the default
% renderer, Painters.  On the other hand, all other renderers make the
% title (etc.) fonts look like shit.  So you had to convert the RGB patches
% in the figure to "color mapping" with rgb2cm.m.
%
% NB: Note that in the titlestring you need to include $s around math
% fonts, but in the xLabels you must not!!

%-------------------------------------------------------------------------%
% Revised: 04/15/14
%   -elimated references to the LaTeX intepreter, b/c the f*ckers at the
%       Mathworks broke it.  The best way to use this now is with
%       matlab2tikz.m
% Created: 10/25/12
%   by JGM
%   -this code was inspired by code written by "us" at matlabcentral
%-------------------------------------------------------------------------%

% ordinary bar plot
h = figure;
bh = bar(data);

% change the bars' colors
ch = get(bh,'children');
cdata = reshape(repmat(barClrs,1,5)',size(barClrs,2),size(barClrs,1)*5)';
cdata = [cdata; ones(1,size(barClrs,2))]; % nan*ones(1,size(barClrs,2))];
set(ch,'facevertexcdata',cdata);
rgb2cm(h);

% modify the maximum value on the plot
ax = axis;
if ~isempty(ymax),  ax(4) = ymax; end
axis(ax);

% label
set(gcf,'Renderer','Painters');
if ~isempty(xlbl), xlabel(xlbl); end
if ~isempty(ylbl), xlabel(ylbl); end
if ~isempty(xTickStrs), set(gca,'XTickLabel',xTickStrs); end
if ~isempty(yTickStrs), set(gca,'YTickLabel',yTickStrs); end
if ~isempty(titlestr), title(titlestr); end


%%% had to kill b/c matlab broke their LaTeX interpreter
% if ~isempty(xlbl), xlabel(xlbl,'Interpreter','Latex','Fontsize',fontsize); end
% if ~isempty(ylbl), ylabel(ylbl,'Interpreter','Latex','Fontsize',fontsize); end
% TexTickLabels(xTickStrs,yTickStrs,fontsize);
% title(titlestr,'Interpreter','Latex','Fontsize',fontsize);


end




