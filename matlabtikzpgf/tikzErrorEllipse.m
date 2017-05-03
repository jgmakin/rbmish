function tikzErrorEllipse(eStats,conf,clrNames,tikzfilename)
% tikzErrorEllipse
% 
% USAGE:
%   tikzErrorEllipse(eStats,conf,clrNames,tikzfilename)
%
% It's much more efficient to have pgf write the ellipses than to have it
% plot a bunch of points that Matlab provided that "happen" to lie on an
% ellipse.
%
% You can change this (in a number of ways), but the biggest ellipse is set
% to be 3in x 3in (??).

%-------------------------------------------------------------------------%
% Revised: 03/26/15
%   -fixed all the scaling issues:
%       --made the largest ellipse take up (approx.) 3x3 in.
%       --placed the coordinate frame just off "'s lower left corner
%       --made the scalebar label set automatically (!)
% Created: 01/21/15
%   by JGM
%-------------------------------------------------------------------------%


%%% TO DO:
% (1) Make the axis labels into arguments!! right now it's all radians rad
%%%%%%%%%%%%%%%%%%%%%%


% r=2 is the number of dimensions (d.o.f.)
confScale = qchisq(conf,2);

% get the "size" (geometric-average radius) of the biggest ellipse
% dd = zeros(size(eStats.Xpct,2),1);
% for iStat = 1:size(eStats.Xpct,2)
%     dd(iStat) = double(det(confScale*eStats.Cvrn(:,:,iStat)));
% end   
% scaleFactor = sqrt(sqrt(max(dd)));
scaleFactor = sqrt(sqrt(max(arrayfun(@(iStat)(...
   det(confScale*gather(eStats.Cvrn(:,:,iStat)))),1:size(eStats.Xpct,2)))));
%%% one sqrt to undo the det, one to move from variance to std. dev. units

% fixed params
axisLength = 0.64;                              % cm
axisLoc = -5/4;                                 % cm
axisLabelLoc = axisLoc + axisLength;            % cm
yrScale = axisLength*scaleFactor;               % radians or whatever


% create the coordinate frame
outtxt = ['\newcommand{\coordinateFrame}[2]{',sprintf('\n\t')];
outtxt = [outtxt,'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- (#1,#2+',num2str(axisLength),');',...
    sprintf('\n\t')];
outtxt = [outtxt,'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- ',sprintf('\n\t\t')];
outtxt = [outtxt,'node[label={[label distance=-40, '];
outtxt = [outtxt,'color=black]{-5}:{',num2str(yrScale,3),'rad}}]{}',...
    sprintf('\n\t\t')];
outtxt = [outtxt,'(#1+',num2str(axisLength),',#2);',sprintf('\n\t')];
outtxt = [outtxt,'}',sprintf('\n')];

% now do for each ellipse...
for iStat = 1:size(eStats.Xpct,2)
    
    Cvrn = eStats.Cvrn(:,:,iStat);
    Xpct = eStats.Xpct(:,iStat);
    M = eStats.N(iStat);
    clr = clrNames{iStat};
    
    % don't plot if the covariance matrix has NaNs in it
    if any(isnan(Cvrn(:)))
        fprintf('covariance has NaN entries --- jgm\n\n');
    else
        
        % plot error ellipse 
        SigmaOneHalf = chol(confScale*Cvrn)'/scaleFactor; % see error_ellipse.m
        Cntr = Xpct/scaleFactor;
        outtxt = [outtxt,...
            tikzEllipse(Cntr,SigmaOneHalf,clr,2.0)];
        
        % plot the standard error (covariance) of the mean, too
        SigmaOneHalf = chol(confScale*Cvrn/M)'/scaleFactor;
        Cntr = Xpct/scaleFactor;
        outtxt = [outtxt,...
            tikzEllipse(Cntr,SigmaOneHalf,clr,2.0)];
    end
end

% now place the coordinate frame 
outtxt = [outtxt,'\coordinateFrame{',num2str(axisLoc),'}'];
outtxt = [outtxt,'{',num2str(axisLoc),'}',sprintf('\n')];
outtxt = [outtxt,'\node[text=black,anchor=west,rotate=90] (yaxislabel) '];
outtxt = [outtxt,'at (',num2str(axisLoc),',',num2str(axisLabelLoc),')'];
outtxt = [outtxt,' {$\prop_2$ (elbow)};',sprintf('\n')];
outtxt = [outtxt,'\node[text=black,anchor=west] (yaxislabel) '];
outtxt = [outtxt,'at (',num2str(axisLabelLoc),',',num2str(axisLoc),')'];
outtxt = [outtxt,' {$\prop_1$ (shoulder)};',sprintf('\n')];

% write the file, wrapped in a tikzpicture
tikzWrite(outtxt,tikzfilename,'[x=2.5cm,y=2.5cm]');

    
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function x = qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.

%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
    dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
    x(s) = x(s)+dx;
    if all(abs(dx) < 1e-6)
        break;
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.

%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
    s = x>0;
    k = 0;
    for jj = 0:n/2-1;
        k = k + (x(s)/2).^jj/factorial(jj);
    end
    F(s) = 1-exp(-x(s)/2).*k;
else
    for ii=1:numel(x)
        if x(ii) > 0
            % F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
            F(ii) = integral(@dchisq,0,x(ii),1e-6,0,n);
        else
            F(ii) = 0;
        end
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));

end
%-------------------------------------------------------------------------%