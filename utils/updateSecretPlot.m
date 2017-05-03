function updateSecretPlot(newdatum,fignum)
% updateSecretPlot  Update a plot with only the current datum
%
% USAGE:
%   updateSecretPlot(newdata,fignum)
%
% If fignum already has data on it, this function just adds one more datum,
% newdata; otherwise, it starts a new plot.

%-------------------------------------------------------------------------%
% Created: 01/22/15
%   by JGM
%-------------------------------------------------------------------------%


set(0,'CurrentFigure',fignum);
dataObjs = get(get(gcf, 'Children'), 'Children');
xdata = get(dataObjs, 'XData');
ydata = get(dataObjs, 'YData');
clf;
if isempty(xdata)
    plot(1,newdatum,'g');
else
    plot([xdata,xdata(end)+1],[ydata,newdatum],'g');
end
%drawnow;

end