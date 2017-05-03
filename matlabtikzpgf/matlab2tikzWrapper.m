function matlab2tikzWrapper(tikzfilename,figHndl,varargin)
% matlab2tikzWrapper    Calls matlab2tikz.m with some standard stuff
%
% USAGE:
%   matlab2tikzWrapper(tikzfilename,figHndl,varargin)
%
% The varargin is for extra colors.

%-------------------------------------------------------------------------%
% Revised: 03/09/17
%   
% Revised: 04/15/14
%   -it's actually better not to tell matlab2tikz a height and width, b/c
%   otherwise it messes up the aspect ratio.  If you want to change the
%   size of the figure, do it in your figures .tex file.  (If you really
%   want to set the aspect ratio to something specific, then go ahead and
%   do it in here.)
% Created: 04/15/14
%   by JGM
%-------------------------------------------------------------------------%




% default dimensions
% defaultDimensionsStr = '\providecommand{\figurewidth}{6cm}\providecommand{\figureheight}{6cm}';

% useful stuff
yrtitle = figHndl.CurrentAxes.Title.String; % extract title (unfortunate)
[pretikzpicture,scaleChanger,xscaleChanger,yscaleChanger,...
    xlabelLocator,ylabelLocator,xticklabelHider,yticklabelHider,...
    legendHider,titleHider] = wrapperbits(yrtitle);

% convert to tikz--with prepended code!!
matlab2tikz([getdir('tikz'),tikzfilename,'.tex'],...
    'figurehandle',figHndl,...
    'parseStrings',false,...
    'extraCode',pretikzpicture,...
    'showInfo', false,...
    'extraAxisOptions',{scaleChanger,xscaleChanger,yscaleChanger,...
    xlabelLocator,ylabelLocator,xticklabelHider,yticklabelHider,...
    legendHider,titleHider},...
    varargin{:}); % ,...
    % 'floatFormat','%.3g');
    % 'showInfo', false,...     
    % 'height', '\figureheight',...
    %'width', '\figurewidth',...
    % 'extraCode',defaultDimensionsStr);

end
