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

% r=2 is the number of dimensions (d.o.f.)
confScale = qchisq(conf,2);

% get the "size" (geometric-average radius) of the biggest ellipse
scaleFactor = sqrt(sqrt(max(arrayfun(@(iStat)(...
    det(confScale*eStats.Cvrn(:,:,iStat))),1:size(eStats.Xpct,2)))));
%%% one sqrt to undo the det, one to move from variance to std. dev. units


% fixed params
axisLength = 0.6;                               % inches
ellipseRadius = 1;                              % inches
axisLoc = -5/4*ellipseRadius;                   % inches
axisLabelLoc = axisLoc + axisLength;            % inches
yrScale = axisLength/ellipseRadius*scaleFactor; % radians or whatever


% text to be written
outtxt = ['\begin{tikzpicture}',carriagereturn];

% create the coordinate frame
outtxt = [outtxt,'\newcommand{\coordinateFrame}[2]{',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- (#1,#2+',num2str(axisLength),'in);',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- ',carriagereturn];
outtxt = [outtxt,sprintf('\t\t'),'node[label={[label distance=-40, '];
outtxt = [outtxt,'color=black]{-5}:{',num2str(yrScale,3),'rad}}]{}',carriagereturn];
outtxt = [outtxt,sprintf('\t\t'),'(#1+',num2str(axisLength),'in,#2);',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'}',carriagereturn];



% now do for each ellipse...
for iStat = 1:size(eStats.Xpct,2)
    
    Cvrn = eStats.Cvrn(:,:,iStat);
    Xpct = eStats.Xpct(:,iStat);
    M = eStats.N(iStat);
    
    if any(isnan(Cvrn(:)))
        fprintf('covariance has NaN entries --- jgm\n\n');
    else
        
        % get transforming matrix via covariance (see error_ellipse.m)
        SigmaOneHalf = chol(confScale*Cvrn)/scaleFactor;
        Xpct = Xpct/scaleFactor;
        
        % make the opacity set-able from the outside w/a command
        outtxt = [outtxt,...
            '\providecommand{\',clrNames{iStat},'opacity}{1}',carriagereturn];
        
        % write into tikz's matrix-scaling thing
        outtxt = [outtxt,'\begin{scope}[cm={',...
            num2str(SigmaOneHalf(1,1)),',',num2str(SigmaOneHalf(1,2)),',',...
            num2str(SigmaOneHalf(2,1)),',',num2str(SigmaOneHalf(2,2)),',(',...
            num2str(Xpct(1)),',',num2str(Xpct(2)),')}]',...
            carriagereturn];
        % outtxt = [outtxt,'\draw [',eStats.clr,',line width=2.0,opacity=\',...
        outtxt = [outtxt,sprintf('\t'),'\draw [',clrNames{iStat},...
            ',line width=2.0,opacity=\',clrNames{iStat},...
            'opacity] (0,0) {} ellipse (',num2str(ellipseRadius),...
            'in and ',num2str(ellipseRadius),'in);',carriagereturn];
        outtxt = [outtxt,'\end{scope}',carriagereturn];
        
        
        
        % Now repeat for the standard error of the mean
        SigmaOneHalf = chol(confScale*Cvrn/M)/scaleFactor;
        Xpct = Xpct/scaleFactor;
        
        % make the opacity set-able from the outside w/a command
        outtxt = [outtxt,...
            '\providecommand{\',clrNames{iStat},'opacity}{1}',carriagereturn];
        
        % write into tikz's matrix-scaling thing
        outtxt = [outtxt,'\begin{scope}[cm={',...
            num2str(SigmaOneHalf(1,1)),',',num2str(SigmaOneHalf(1,2)),',',...
            num2str(SigmaOneHalf(2,1)),',',num2str(SigmaOneHalf(2,2)),',(',...
            num2str(Xpct(1)),',',num2str(Xpct(2)),')}]',...
            carriagereturn];
        % outtxt = [outtxt,'\draw [',eStats.clr,',line width=2.0,opacity=\',...
        outtxt = [outtxt,sprintf('\t'),'\draw [',clrNames{iStat},...
            ',line width=2.0,opacity=\',clrNames{iStat},...
            'opacity] (0,0) {} ellipse (',num2str(ellipseRadius),...
            'in and ',num2str(ellipseRadius),'in);',carriagereturn];
        outtxt = [outtxt,'\end{scope}',carriagereturn];
        
        
    end
end


% now place the coordinate frame 
outtxt = [outtxt,'\coordinateFrame{',num2str(axisLoc),'in}'];
outtxt = [outtxt,'{',num2str(axisLoc),'in}',carriagereturn];
outtxt = [outtxt,'\node[text=black,anchor=west,rotate=90] (yaxislabel) '];
outtxt = [outtxt,'at (',num2str(axisLoc),'in,',num2str(axisLabelLoc),'in)'];
outtxt = [outtxt,' {$\prop_2$ (elbow)};',carriagereturn];
outtxt = [outtxt,'\node[text=black,anchor=west] (yaxislabel) '];
outtxt = [outtxt,'at (',num2str(axisLabelLoc),'in,',num2str(axisLoc),'in)'];
outtxt = [outtxt,' {$\prop_1$ (shoulder)};',carriagereturn];

% close off the picture
outtxt = [outtxt,'\end{tikzpicture}%'];

% write to file
fid = fopen(outfile,'wt+');
fwrite(fid,outtxt);
fclose(fid);

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
            F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
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