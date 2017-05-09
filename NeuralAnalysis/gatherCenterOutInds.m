function cntrOutInds = gatherCenterOutInds(...
    cntrdists,distThr,mvmtspeeds,speedThr,T,TOPLOT,X,xmin,xmax)
% Extract *just* the center-out portions of a reach trajectory.  A speed
% threshold (speedThr) is used to determine when a movement has started or
% ended.  A distance threshold (distThr) is used to make sure the movement
% finishes "far" from the center.  (NB: This will also throw out out-center
% reaches, which is what you want.)  Hence the function also needs to know
% the distance from the center (cntrdists) and the speeds (mvmtspeeds) at
% all point along the reach.  
%
% The function also takes an argument T, in case one wants to pad the
% reaches (of length size(cntrdists,1)) with zeros at the end.  This is
% useful for putting the result into a matrix/tensor, rather than a cell
% array.
%
% The last three arguments, X, xmin, and xmax, are for plotting purposes
% only, and are only invoked if TOPLOT==1.


%-------------------------------------------------------------------------%
% Cribbed: 04/08/15
%   from extractCenterOutInds
%   by JGM
%-------------------------------------------------------------------------%


% find where the center-out reach starts and ends
belowThreshTimes = zeros(1,T);
belowThreshTimes(mvmtspeeds<speedThr) = 1;
belowThreshChanges = diff(belowThreshTimes);

% get onset and offset indices for all movements
%%% FUNCTIONIZE ME
onsets = [belowThreshChanges 0]==-1;
onsetInds = find(onsets);
offsets = [0 belowThreshChanges]==1;
offsetInds = find(offsets);
Nmvmts = length(onsetInds);
if Nmvmts>length(offsetInds)
    offsetInds = [offsetInds, T];
    offsets(T) = 1;
end

% make sure the final position of these movements is away from the cntr
cntrOutOnsetInds = onsetInds(cntrdists(offsetInds) > distThr);
cntrOutOffsetInds = offsetInds(cntrdists(offsetInds) > distThr);
%%% you could also have required that the *start* of each trajectories
%%% be *close* to the center


% plot
if TOPLOT
    figure(5); clf; hold on;
    plot(mvmtspeeds)
    scatter(onsetInds,mvmtspeeds(onsetInds),'g')
    scatter(offsetInds,mvmtspeeds(offsetInds),'r')
    scatter(cntrOutOnsetInds,mvmtspeeds(cntrOutOnsetInds),'g*')
    scatter(cntrOutOffsetInds,mvmtspeeds(cntrOutOffsetInds),'r*')
    hold off;
    pause()
    
    for iMvmt = 1:length(cntrOutOnsetInds)
        figure(6); clf; hold on;
        xx = squeeze(X(1,cntrOutOnsetInds(iMvmt):cntrOutOffsetInds(iMvmt)));
        yy = squeeze(X(2,cntrOutOnsetInds(iMvmt):cntrOutOffsetInds(iMvmt)));
        plot(xx,yy);
        scatter(xx(end),yy(end),'rx');
        axis([xmin(1),xmax(1),xmin(2),xmax(2)]);
        hold off;
        pause()
    end
end


% store
cntrOutInds = [cntrOutOnsetInds;cntrOutOffsetInds];

end
%-------------------------------------------------------------------------%