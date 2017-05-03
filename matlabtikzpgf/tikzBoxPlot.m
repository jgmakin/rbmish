function tikzBoxPlot(x,y,clrNames,xticklabels,ylabel,optaxisargs,...
    optplots,tikzfilename)
% tikzBoxPlot   Write out a multiple-boxplot tikz file

%-------------------------------------------------------------------------%
% Created: 09/18/15
%   by JGM
%-------------------------------------------------------------------------%


% params
Ngroups = length(y);
Noptaxisargs = length(optaxisargs);
Noptplots = length(optplots);


% ----------- init
outtxt = ['\usepgfplotslibrary{statistics}%',sprintf('\n')];
outtxt = [outtxt,'\input{\texdir/ucsf/REFH/REFH.sty}%',sprintf('\n')];
outtxt = [outtxt,'\providecommand{\thisTikzPicScale}{1.0}%',sprintf('\n\t')];
% ----------- axis options
outtxt = [outtxt,'\begin{axis}',sprintf('\n\t\t')];
outtxt = [outtxt,'[',sprintf('\n\t\t')];
outtxt = [outtxt,'boxplot/draw direction=y,',sprintf('\n\t\t')];
outtxt = [outtxt,'xtick={',num2str(x(1)),','];
outtxt = [outtxt,num2str(x(2)),',...,',num2str(x(end)),'},',sprintf('\n\t\t')];
outtxt = [outtxt,'xticklabels={'];
for iTick = 1:length(xticklabels)
    outtxt = [outtxt,xticklabels{iTick},','];
end
outtxt = [outtxt,'},',sprintf('\n\t\t')];
outtxt = [outtxt,'clip=false,',sprintf('\n\t\t')];
outtxt = [outtxt,'ylabel={',ylabel,'},',sprintf('\n\t\t')];
outtxt = [outtxt,'scale=\thisTikzPicScale,',sprintf('\n\t\t')];
outtxt = [outtxt,'every axis x label/.style={at={(ticklabel cs:0.5)},anchor=north},',sprintf('\n\t\t')];
outtxt = [outtxt,'every axis y label/.style={at={(ticklabel cs:0.5)},anchor=south,rotate=90},',sprintf('\n\t\t')];
% ------------ optional axis options
for iOpt = 1:Noptaxisargs
    outtxt = [outtxt,optaxisargs{iOpt},',',sprintf('\n\t\t')];
end
outtxt = [outtxt,']',sprintf('\n\t\t')];
%------------ boxplot 
for iGroup = 1:Ngroups
    thisY = y{iGroup};
    outtxt = [outtxt,'\addplot['];
    % ------- plot properties
    outtxt = [outtxt,'color=',clrNames{iGroup},','];
    outtxt = [outtxt,'line width=1.0,solid,',];
    outtxt = [outtxt,'forget plot,',];
    outtxt = [outtxt,'mark=*,',];
    outtxt = [outtxt,'boxplot={draw position=',num2str(x(iGroup)),'}'];
    outtxt = [outtxt,']',sprintf('\n\t\t\t')];
    outtxt = [outtxt,'table[row sep=\\, y index=0]{%',sprintf('\n\t\t\t')];
    outtxt = [outtxt,'data\\',sprintf('\n')];
    for iDatum = 1:length(thisY)
        outtxt = [outtxt,num2str(thisY(iDatum)),'\\ '];
        %%% might want to tell num2str how many significant digits?
    end
    outtxt = [outtxt,sprintf('\n\t\t'),'};',sprintf('\n')];
end 
%------------ extra stuff
for iOpt = 1:Noptplots
    outtxt = [outtxt,optplots{iOpt},sprintf('\n')];
end
outtxt = [outtxt,sprintf('\t'),'\end{axis}',sprintf('\n')];

% write out as tikz file, wrapped in a tikzpicture
tikzWrite(outtxt,tikzfilename);


end