function ax = getCustomAxesPos(rows,cols,space)
%GETCUSTOMAXESPOS   Creates a matrix of handles to custom subplot axes
%   AX = GETCUSTOMAXESPOS(ROWS,COLS,SPACE) creates a [ROWS x COLS] matrix
%   AX of graphics handles to subplots with SPACE percent gray space
%   between them.

% Finds the dimensions of the figure as a whole by creating one big axes
% (which we will delete later) when we do not need it any longer.
BigAx = newplot;
set(BigAx,'Visible','off','color','none')
pos = get(BigAx,'Position');

% Calcuates the width and height of each subplot.
% Also where you can set the percent spacing between axes.
width = pos(3)/cols;            % width of each subplot
height = pos(4)/rows;           % height of each subplot
% width = height;
%space = 0;                      % percent space between axes (.02 = 2%)
pos(1:2) = pos(1:2) + space*[width height];

% Deleting the big axes
BigAxParent = get(BigAx,'Parent');
delete(BigAx);

% Makes handles to each of the subplot axes. Sets the position of each
% subplot to eliminate all the gray spaces between them. The handles to
% each subplot can be found in matrix 'ax'.
ax = zeros(rows,cols);
for i=rows:-1:1,
    for j=cols:-1:1,
        axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height ...
            width*(1-space) height*(1-space)];
        ax(i,j) = axes('Position',axPos,'parent',BigAxParent);
        % colormap('gray');
        % set(ax(i,j),'visible','on'); % You can use this to make the axes
        % visible, we commented this out for now.
        
    end
end