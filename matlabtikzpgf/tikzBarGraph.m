function tikzBarGraph(xdata,ydata,errorBars,ymin,xticklabels,...
    xlabel,ylabel,title,clrNames,height,barwidth,plottype,legendentries,...
    tikzfilename,varargin)
% This function was unfortunately necessitated by the incompatibility of
% your workarounds to get bar plots with a different color on every bar (on
% the one hand) and matlab2tikz (on the other).

%-------------------------------------------------------------------------%
% Revised: 04/27/17
%   -finally got a satisfactory setting of widths for group plots (by
%   taking the bar spacing into account)
%   -altered to allow for extra axis options as well as extra code before
%   the plots are closed (before \end{axis})
% Revised: 06/06/16
%   -fixed a bug that left "xticklabels=" rather than nothing at all when
%   the input xticklabels is blank.  Ditto for "xtick=".
% Revised: 03/16/16
%   -added new first argument, xdata, for easier use with numerically
%   labelled bars
%   -allowed the xticklabels to be empty, in which case the function lets
%   pgf take care of the xticks and xticklabels
% Revised: 11/04/15
%   -changed to allow multiple bar graphs on top of each other; now ydata
%   etc. can be matrices rather than vectors
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
extraAxisText = defaulter('extraAxisText',[],varargin{:});
extraPlotText = defaulter('extraPlotText',[],varargin{:});
extraTikzText = defaulter('extraTikzText',[],varargin{:});


% these are useful bits of code, in string form
[pretikzpicture,scaleChanger,xscaleChanger,yscaleChanger,...
    xlabelLocator,ylabelLocator,xticklabelHider,yticklabelHider,...
    legendHider,titleHider] = wrapperbits(title);

% adjust the style definition for y tick labels to set the precision
yticklabelHider = sprintf([...
    '%s\t',...
    '/pgf/number format/fixed,%%\n\t',...
    '/pgf/number format/precision=5,%%\n',...
    '}'],yticklabelHider(1:end-1));
%%% evidently you don't want this done in matlab2tikzWrapper, which also
%%% calls wrapperbits.m, otherwise you'd place it in there


% figure params
[Ny,Ngrps] = size(ydata);
incr = min(diff(xdata))/2;
xmin = xdata(1) - incr;
xmax = xdata(end) + incr;
ymax = 1.1*max(ydata(:));
widthscale = 1.3*length(xdata)*Ngrps;
%%% by default, which you accept, the spacing between bars is 2pt=2/72 in.


% build x ticks, x-tick labels, and colors
if isempty(xticklabels)
    xtickStr = '';
    xticklabelsStr = '';
else
    xtickStr = 'xtick={';
    xticklabelsStr = 'xticklabels={';
    for xind = 1:length(xdata)
        xtickStr = sprintf('%s%g,',xtickStr,xdata(xind));
        xticklabelsStr = sprintf('%s{%s},',xticklabelsStr,xticklabels{xind});
        %     thisColor = clrs(xind,:);
        %     defineCustomColorsStr = [defineCustomColorsStr,...
        %         '\definecolor{mycolor',...
        %         num2str(xind),...
        %         '}{rgb}{',...
        %         num2str(thisColor(1)),',',...
        %         num2str(thisColor(2)),',',...
        %         num2str(thisColor(3)),'}',...
        %         sprintf('\n')];
    end
    % get rid of last comma
    xtickStr = sprintf('%s},\n',xtickStr(1:end-1));
    xticklabelsStr = sprintf('%s},\n',xticklabelsStr(1:end-1));
end



% enclose axis labels and title in curly braces
xlabel = ['{',xlabel,'}'];
ylabel = ['{',ylabel,'}'];

% text to be written
outtxt = [pretikzpicture,newline,...
    sprintf('\\providecommand{\\yrbarwidth}{\\thisTikzPicScale*%g in}%%\n',barwidth),...
    sprintf('\\begin{axis}[%%\n'),...
    sprintf('height=%g in,\n',height),...
    sprintf('area legend,\n'),...
    sprintf('scale only axis,\n'),...
    sprintf('xlabel=%s,\n',xlabel),...
    sprintf('ymin=%d,\n',ymin),...
    sprintf('ymax=%d,\n',ymax),...
    sprintf('ylabel=%s,\n',ylabel),...
    sprintf('width=\\thisTikzPicScale*%g*(\\yrbarwidth + 2pt),\n',widthscale)];
%%% Ideally, the 2pt would be replaced by a reference to ybar....
%%%outtxt = [outtxt,sprintf(['axis x line*=bottom,\n','axis y line*=left,\n'])];
outtxt = [outtxt,sprintf('axis x line*=bottom,\naxis y line*=left,\n')];
outtxt = [outtxt,sprintf('x axis line style={draw=none},\n')];
outtxt = [outtxt,sprintf('extra y ticks= 0,\nextra y tick labels=,\n')];
outtxt = [outtxt,sprintf('extra y tick style= {grid=major},\n')];
outtxt = [outtxt,sprintf('%s,\n%s,\n%s,\n%s,\n%s,\n%s,\n%s,\n%s,\n%s,\n',...
    scaleChanger,...
    xscaleChanger,...
    yscaleChanger,...
    xlabelLocator,...
    ylabelLocator,...
    xticklabelHider,...
    yticklabelHider,...
    legendHider,...
    titleHider)];

switch plottype
    case 'grouped'
        
        % finish axis definition
        outtxt = [outtxt,...
            sprintf('xtick = data,\n'),...
            sprintf('symbolic x coords=%s',xticklabelsStr(13:end)),...
            sprintf('ybar,\n'),...
            sprintf('bar width=\\yrbarwidth,\n'),...
            sprintf('major x tick style = transparent,\n'),...
            sprintf('enlarge x limits={abs={%g*(\\yrbarwidth + 2pt)}},\n',Ngrps/2),...
            extraAxisText,...
            sprintf(']\n\n')];
        
        % loop across groups...
        for iGrp = 1:Ngrps
            outtxt = [outtxt,...
                plotOneBarGroup(xticklabels,ydata(:,iGrp),...
                clrNames{1,iGrp},errorBars(:,:,iGrp))];
        end
        
        % legend
        if ~isempty(legendentries)
            legendStr = '\legend{';
            for iGrp = 1:Ngrps
                legendStr = sprintf('%s%s, ',legendStr,legendentries{iGrp});
            end
            legendStr = sprintf('%s}',legendStr);
            
            outtxt = [outtxt,sprintf('%s\n',legendStr)];
        end
        
        
    case 'overlapping'
        
        
        % finish axis definition
        outtxt = [outtxt,...
            sprintf('xmin=%d,\n\t',xmin),...
            sprintf('xmax=%d,\n\t',xmax),...
            sprintf('%s%s',xtickStr,xticklabelsStr),...
            extraAxisText,...
            sprintf(']\n\n')];
        
        % make plots
        if Ngrps > 1, opacity = 0.5; else opacity = 1; end
        for iGrp = 1:Ngrps
            for xind = 1:Ny
                outtxt = [outtxt,...
                    plotOneBar(xdata(xind),ydata(xind,iGrp),...
                    sprintf('%s,opacity=%.2f',clrNames{xind,iGrp},opacity),...
                    errorBars(xind,:,iGrp))];
            end
        end
        
        
    otherwise
        fprintf('defaulting to overlapping bars');
        
        % finish axis definition
        outtxt = [outtxt,...
            sprintf('xmin=%d,\n\t',xmin),...
            sprintf('xmax=%d,\n\t',xmax),...
            sprintf('%s%s',xtickStr,xticklabelsStr),...
            extraAxisText,...
            sprintf(']\n\n')];
        
        % make plots
        if Ngrps > 1, opacity = 0.5; else opacity = 1; end
        for iGrp = 1:Ngrps
            for xind = 1:Ny
                outtxt = [outtxt,...
                    plotOneBar(xdata(xind),ydata(xind,iGrp),...
                    sprintf('%s,opacity=%.2f',clrNames{xind,iGrp},opacity),...
                    errorBars(xind,:,iGrp))];
            end
        end
end



% the end
outtxt = [outtxt,newline,extraPlotText,newline,'\end{axis}'];

% write to file as tikzpicture
tikzWrite(outtxt,tikzfilename,extraTikzText)

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
function outtxt = plotOneBar(xdatum,ydatum,clrName,errorBar)

outtxt = [...
    sprintf('\\addplot[ybar,bar width=\\yrbarwidth,draw=none,fill=%s]\n',clrName),...
    sprintf('\tplot[error bars/.cd, y dir=both, y explicit]\n'),...
    sprintf('\ttable[row sep=crcr')];
if any(errorBar)
    outtxt = [outtxt,...
        ', y error plus index=2, y error minus index = 3] {',...
        sprintf('%d %0.10d %0.10d %0.10d',xdatum,ydatum,errorBar(1),errorBar(2))];
else
    outtxt = [outtxt,sprintf('] {%d %0.10d',xdatum,ydatum)];
end
outtxt = [outtxt,sprintf('\\};\n')];
% 'mycolor',num2str(xind),...
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function outtxt = plotOneBarGroup(xticklabels,ydata,clrName,errorBars)

% use error bars?
USEERRORBARS = ~any(isnan(errorBars(:)));

% \addplot...
outtxt = sprintf('\\addplot[draw=none,fill=%s]\n',clrName);
if USEERRORBARS
    outtxt = [outtxt,sprintf('\tplot[error bars/.cd, y dir=both, y explicit]\n')];
end
outtxt = [outtxt,sprintf('\tcoordinates {\n')];

% content of the plot (doesn't hurt to write the error bars)
for xind = 1:length(xticklabels)
    outtxt = [outtxt,...
        sprintf('\t\t(%s, %d)',xticklabels{xind},ydata(xind)),...
        sprintf('-= (0,%d)',-errorBars(xind,1)),...
        sprintf('+= (0,%d)', errorBars(xind,2)),...
        sprintf('\n')];
end
outtxt = [outtxt,sprintf('\t};\n')];

end
%-------------------------------------------------------------------------%































