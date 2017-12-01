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



 %%%xticklabel style={/pgf/number format/frac}
 
% default dimensions
% defaultDimensionsStr = '\providecommand{\figurewidth}{6cm}\providecommand{\figureheight}{6cm}';

% useful stuff
yrtitle = figHndl.CurrentAxes.Title.String; % extract title (unfortunate)
[defaultExtraCode,defaultExtraAxisOptions] = wrapperbits(yrtitle);

% if there is *extra* extra code, etc. for matlab2tikz.m, pull it out
variable_args = varargin;
[extraCode, variable_args] = defaulterplus('extraCode',...
    defaultExtraCode,variable_args{:});
[extraAxisOptions, variable_args] = defaulterplus('extraAxisOptions',...
    defaultExtraAxisOptions,variable_args{:});

% convert to tikz--with prepended code!!
matlab2tikz([getdir('tikz'),tikzfilename,'.tex'],...
    'figurehandle',figHndl,...
    'parseStrings',false,...
    'extraCode',extraCode,...
    'showInfo', false,...
    'extraAxisOptions',extraAxisOptions,...
    variable_args{:}); % ,...
    % 'floatFormat','%.3g');
    % 'showInfo', false,...     
    % 'height', '\figureheight',...
    %'width', '\figurewidth',...
    % 'extraCode',defaultDimensionsStr);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [val_cell, variable_args] = defaulterplus(name,val_cell,varargin)
% instead of replacing with the default, tack it on to the default

iVal = find(strcmp(varargin,name));
if ~isempty(iVal)
    new_vals = varargin{iVal+1};
    val_cell = [val_cell, new_vals];
    variable_args = varargin([1:iVal-1,iVal+2:end]);
else
    variable_args = varargin;
end

end
%-------------------------------------------------------------------------%


















