function [pretikzpicture,scaleChanger,xscaleChanger,yscaleChanger,...
    xlabelLocator,ylabelLocator,xticklabelHider,yticklabelHider,...
    legendHider,titleHider] = wrapperbits(yrtitle)
% wrapperbits   Strings that will be written at the beginning of tikz file

%-------------------------------------------------------------------------%
% Created: 03/29/17
%   by JGM
%-------------------------------------------------------------------------%


pretikzpicture = sprintf([...
    '\\providecommand{\\thisTikzPicScale}{1}%%\n',...
    '\\providecommand{\\thisXscale}{1}%%\n',...
    '\\providecommand{\\thisYscale}{1}%%\n',...
    '\\providecommand{\\thisXlabelopacity}{1}%%\n',...
    '\\providecommand{\\thisYlabelopacity}{1}%%\n',...
    '\\providecommand{\\thisXticklabelopacity}{1}%%\n',...
    '\\providecommand{\\thisYticklabelopacity}{1}%%\n',...
    '\\provideboolean{SKIPLEGEND}%%\n',...
    '\\provideboolean{SKIPTITLE}%%',... % m2t takes care of last \n
    ]);


% extra axis options
scaleChanger = 'scale=\thisTikzPicScale';
xscaleChanger = 'x post scale=\thisXscale';
yscaleChanger = 'y post scale=\thisYscale';
xlabelLocator = sprintf([...
    'every axis x label/.style={%%\n\t',...
    'at={(xticklabel cs:0.5)},%%\n\t',...
    'anchor=north,%%\n\t',...
    'opacity=\\thisXlabelopacity,%%\n',...
    '}']);
%%% 'at={($(yticklabel cs:0) + (xticklabel cs:0.5)$)},%%\n\t',...
ylabelLocator = sprintf([...
    'every axis y label/.style={%%\n',...
    '\tat={(ticklabel cs:0.5)},%%\n',...
    '\tanchor=south,%%\n',...
    '\trotate=90,%%\n',...
    '\topacity=\\thisYlabelopacity,%%\n',...
    '}']);
xticklabelHider = sprintf([...
    'every x tick label/.style={%%\n',...
    '\topacity=\\thisXticklabelopacity,%%\n',...
    '}']);
yticklabelHider = sprintf([...
    'every y tick label/.style={%%\n',...
    '\topacity=\\thisYticklabelopacity,%%\n',...
    '}']);
legendHider = sprintf([...
    'every axis legend/.append code={%%\n\t',...
    '\\ifthenelse{\\boolean{SKIPLEGEND}}{%%\n\t\t'...
    '\\renewcommand\\addlegendentry[2][]{}%%\n\t',...
    '}{}%%\n}']);
titleHider = sprintf([...
    'title={%%\n\t',...
	'\\ifthenelse{\\boolean{SKIPTITLE}}{%%\n\t',...
    '}{%%\n\t\t%s%%\n\t}%%\n}',...
    ],yrtitle);

end