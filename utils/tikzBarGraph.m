function tikzBarGraph(ydata,errorBars,height,barwidth,xticklabels,xlabel,...
    ylabel,title,clrNames,tikzfilename)
% This function was unfortunately necessitated by the incompatibility of
% your workarounds to get bar plots with a different color on every bar (on
% the one hand) and matlab2tikz (on the other).

%-------------------------------------------------------------------------%
% Revised: 01/23/14
%   -changed to use error bars (which required a new second argument)
%   -changed
%       outfile = [outfile,fafasd]
%       outfile = [outfile,gafd]
%   to 
%       outfile = [outfile,fafasd,...
%           gafd,...
%   etc.
% Revised: 10/22/14
%   -changed matrix of colors to cell array of color names (these should be
%   defined in rbmish.sty---or etc.)
% Created: 04/16/14
%   by JGM
%-------------------------------------------------------------------------%


%%% TO DO
%%% (1) height??
%%% (2) fix the height to accommodate the error bars!!
%%% (3) do something, here or elsewhere, to make no error bars if errorBars
%%% is all zeros



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



% you will need the carriagereturn 
load ../toys/filez/carriage.mat
outfile = [yrtikzdir,tikzfilename,'.tex'];


% figure params
xmax = length(ydata)+0.5;
ymin = 0;
ymax = max(ydata) + 0.1*max(ydata); % 0.0016;
% barwidth = 0.32;  % inches
width = barwidth*(xmax+1);  % 4.52;
% height = 3.57;


% build x ticks, x-tick labels, and colors
xtickStr = '{';
xticklabelsStr = '{';
for xind = 1:length(ydata)
    xtickStr = [xtickStr,num2str(xind),',']; 
    xticklabelsStr = [xticklabelsStr,'{',xticklabels{xind},'},'];

    
%     thisColor = clrs(xind,:);
%     defineCustomColorsStr = [defineCustomColorsStr,...
%         '\definecolor{mycolor',...
%         num2str(xind),...
%         '}{rgb}{',...
%         num2str(thisColor(1)),',',...
%         num2str(thisColor(2)),',',...
%         num2str(thisColor(3)),'}',...
%         carriagereturn];
end
xtickStr = [xtickStr(1:end-1),'}'];                 % get rid of the last
xticklabelsStr = [xticklabelsStr(1:end-1),'}'];     % comma


% enclose axis labels and title in curly braces
xlabel = ['{',xlabel,'}'];
ylabel = ['{',ylabel,'}'];
title = ['{',title,'}'];



% text to be written
outtxt = ['\begin{tikzpicture}',carriagereturn,...
    '\providecommand{\thisTikzPicScale}{1}%',carriagereturn,...
    '\begin{axis}[%',carriagereturn,...
    sprintf('\t'),'width=\thisTikzPicScale*',num2str(width),'in,',carriagereturn,...
    sprintf('\t'),'height=',num2str(height),'in,',carriagereturn,...
    sprintf('\t'),'area legend,',carriagereturn,...
    sprintf('\t'),'scale only axis,',carriagereturn,...
    sprintf('\t'),'xmin=0.5,',carriagereturn,...
    sprintf('\t'),'xmax=',num2str(xmax),',',carriagereturn,...
    sprintf('\t'),'xtick=',xtickStr,',',carriagereturn,...
    sprintf('\t'),'xticklabels=',xticklabelsStr,',',carriagereturn,...
	sprintf('\t'),'xlabel=',xlabel,',',carriagereturn,...
    sprintf('\t'),'ymin=',num2str(ymin),',',carriagereturn,...
    sprintf('\t'),'ymax=',num2str(ymax),',',carriagereturn,...
    sprintf('\t'),'ylabel=',ylabel,',',carriagereturn,...
    sprintf('\t'),'title=',title,',',carriagereturn,...
    sprintf('\t'),'axis x line*=bottom,',carriagereturn,...
    sprintf('\t'),'axis y line*=left,',carriagereturn,...
    sprintf('\t'),'scale=\thisTikzPicScale,',carriagereturn,...
    sprintf('\t'),'ylabel absolute, ylabel style={yshift=-1em}',carriagereturn,...
    ']',carriagereturn,carriagereturn];

% the bars themselves: first and last (for some reason)...
% outtxt = [outtxt...
%     '\addplot[color=black,solid,forget plot] table[row sep=crcr] {0 0\\'...
%     num2str(xmax),...
%     '0\\};',carriagereturn];

% ...the rest
for xind = 1:length(ydata)
    outtxt = [outtxt,'\addplot[ybar,bar width=',num2str(barwidth),...
        'in,draw=black,fill=',clrNames{xind},']',carriagereturn,...
        sprintf('\t'),'plot[error bars/.cd, y dir=both, y explicit]',...
        carriagereturn,sprintf('\t'),'table[row sep=crcr'];
    if errorBars(xind,1)||errorBars(xind,2)
        outtxt = [outtxt,...
            ', y error plus index=2, y error minus index = 3] {',...
            num2str(xind),' ',num2str(ydata(xind),10),' '...
            num2str(errorBars(xind,1),10),' ',num2str(errorBars(xind,2),10)];
    else
        outtxt = [outtxt,'] {',num2str(xind),' ',num2str(ydata(xind),10)];
    end
    outtxt = [outtxt,'\\};',carriagereturn];
    % 'mycolor',num2str(xind),...
end

% the end 
outtxt = [outtxt,carriagereturn,'\end{axis}',carriagereturn,...
    '\end{tikzpicture}%'];


% write to file
fid = fopen(outfile,'wt+');
fwrite(fid,outtxt);
fclose(fid);

end