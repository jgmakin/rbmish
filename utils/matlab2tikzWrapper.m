function matlab2tikzWrapper(tikzfilename,figHndl,varargin)
% You will generally use the amazing function matlab2tikz.m with these bits
% of code first.

%-------------------------------------------------------------------------%
% Revised: 04/15/14
%   -it's actually better not to tell matlab2tikz a height and width, b/c
%   otherwise it messes up the aspect ratio.  If you want to change the
%   size of the figure, do it in your figures .tex file.  (If you really
%   want to set the aspect ratio to something specific, then go ahead and
%   do it in here.)
% Created: 04/15/14
%   by JGM
%-------------------------------------------------------------------------%

% directory
[blank, name] = system('hostname');
switch strtrim(name)
    case 'kobayashi-maru'
        yrtikzdir = 'C:\Documents and Settings\makin\My Documents\#texs\tikzpics\';
    case {'CUPCAKE','Themistocles'}
        yrtikzdir = 'C:\Users\makin\Documents\#texs\tikzpics\';
    case {'mushroom','keck-phaser1','pepperoni','zamfir'}
        yrtikzdir = 'C:\Users\makin\Documents\';
    case 'domestica'
        yrtikzdir = '~/tikzpics/';
  otherwise
        error('unknown host -- jgm');
end

% default dimensions
% defaultDimensionsStr = '\providecommand{\figurewidth}{6cm}\providecommand{\figureheight}{6cm}';
defaultScaleStr = '\providecommand{\thisTikzPicScale}{1}';

% colors....

% convert to tikz--with prepended code!!
matlab2tikz([yrtikzdir,tikzfilename,'.tex'],...
    'figurehandle',figHndl,...
    'parseStrings',false,...
    'extraAxisOptions','scale=\thisTikzPicScale',...
    'extraCode',defaultScaleStr,...
    'showInfo', false,...
    'extraColors',varargin);
    % 'showInfo', false,...     
    % 'height', '\figureheight',...
    %'width', '\figurewidth',...
    % 'extraCode',defaultDimensionsStr);

end
