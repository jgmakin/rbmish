function outtxt = tikzEllipse(Cntr,SigmaOneHalf,clr,linewidth)
% tikzEllipse   String of tikz code for an ellipse
%
% USAGE:
%   outtxt = tikzEllipse(Cntr,SigmaOneHalf,clr,linewidth)
%
% Given an ellipse described by a square-root matrix, SigmaOneHalf, and 
% centered at Cntr, provides tikzcode to plot it with color clr and width
% linewidth.

%-------------------------------------------------------------------------%
% Revised: 02/29/16
%   -removed location and ellipse radius arguments, which were really
%   superfluous
% Created: 02/27/16
%   -by JGM
%-------------------------------------------------------------------------%

% make the opacity set-able from the outside w/a command
outtxt = ['\providecommand{\',clr,'opacity}{1}',sprintf('\n')];

% write into tikz's matrix-scaling thing
outtxt = [outtxt,'\begin{scope}[cm={',...
    num2str(SigmaOneHalf(1,1)),',',num2str(SigmaOneHalf(1,2)),',',...
    num2str(SigmaOneHalf(2,1)),',',num2str(SigmaOneHalf(2,2)),',(',...
    num2str(Cntr(1)),',',num2str(Cntr(2)),')}]',...
    sprintf('\n')];
outtxt = [outtxt,sprintf('\t'),'\draw [',clr,...
    ',line width=',num2str(linewidth),',opacity=\',clr,...
    'opacity] (0,0) {} ellipse (1 and 1);',sprintf('\n')];
outtxt = [outtxt,'\end{scope}',sprintf('\n')];

end