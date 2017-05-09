function TexTickLabels(xTickStrs,yTickStrs,fontsize)

%-------------------------------------------------------------------------%
% Revised: 11/02/12
%   -restored ability to modify ylabels
%   -modified so that empty args *TickStr tells matlab to use the original
%   tick labels (but in LaTeX font, in text boxes)
% Revised: 10/25/12
%   -made into a function
% Revised: 05/20/11
%   -made it work
% Cribbed: 05/20/11
%   from the Mathworks' website
%   by JGM
%-------------------------------------------------------------------------%

% get the spacing right
set(gca,'FontSize',fontsize,'Fontname','Times-Roman')

% Get tick mark positions
yTicks = get(gca,'ytick');
xTicks = get(gca,'xtick');

% get length of y axis
yLength = yTicks(end) - yTicks(1);

% get some more useful params
ax = axis;                                      % get leftmost x-position
HorizontalOffset = 0.1;
VerticalOffset = -yLength/20;

% Reset the xtick labels in desired font
if ~isempty(xTickStrs)
    set(gca,'xticklabel',[])                    % remove tick labels
    for xx = 1:length(xTicks)
        
        % get tick label
        if length(xTickStrs) < xx,
            TickLabel = num2str(xTicks(xx));
        else
            TickLabel = xTickStrs{xx};
        end
        
        % Create text box and set appropriate properties
        MakeTextBox(xTicks(xx),yTicks(1)+VerticalOffset,TickLabel,'Center',fontsize);
    end
end


% Reset the ytick labels in desired font
if ~isempty(yTickStrs)
    set(gca,'yticklabel',[])                    % Remove tick labels
    for yy = 1:length(yTicks)
        
        % get tick label
        if length(yTickStrs) < yy,
            TickLabel = num2str(yTicks(yy));
        else
            TickLabel = yTickStrs{yy};
        end
        
        % Create text box and set appropriate properties
        MakeTextBox(ax(1)-HorizontalOffset,yTicks(yy),TickLabel,'Right',fontsize);
    end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function MakeTextBox(xLOC,yLOC,str,halign,fontsize)

% h = zeros(length(xTickStrs));

text(xLOC,yLOC,['$', str, '$'],'HorizontalAlignment',halign,...
    'Interpreter', 'latex','FontSize',fontsize,'FontWeight','bold',...
    'FontName','FixedWidth');

end
%-------------------------------------------------------------------------%