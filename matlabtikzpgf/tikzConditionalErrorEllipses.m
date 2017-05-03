function tikzConditionalErrorEllipses(Nx,Ny,scale,wts,params)
% tikzConditionalErrorEllipses  Plot conditional error covariances in tikz
% 
% USAGE:
%   tikzConditionalErrorEllipses(15,15,2.5,wts,params);
% See e.g. MIvDE.m.
%
% Reproduces Fig. 2a from Makin2013b, "Learning Multisensory Integration
% and Coordinate Transformation via Density Estimation"

%-------------------------------------------------------------------------%
% Created: 02/29/16
%   -by JGM
%-------------------------------------------------------------------------%


% r=2 is the number of dimensions (d.o.f.)
conf = 0.95;
confScale = qchisq(conf,2);

% init
Nexamples = 8000;
smin = params.smin(:,strcmp(params.mods,params.NS));
smax = params.smax(:,strcmp(params.mods,params.NS));
S = grid2world([1 1; params.N ,params.N],[smin,smax],params);


% points at which the decoder will be evaluated
s1 = linspace(smin(1),smax(1),Nx);
s2 = linspace(smin(2),smax(2),Ny);


% make a node in the lower corner
outtxt = ['\node[anchor=south west] (nodeOne) at ('...
    num2str(S(1,1)),',',num2str(S(1,2)),') {};',sprintf('\n')];

% create an axis
outtxt = [outtxt,'\begin{axis}[',sprintf('\n\t'),...
	'at={(nodeOne.south west)},',sprintf('\n\t'),...
	'scale only axis,',sprintf('\n\t'),...
    'xtick={-1.5707963267949,-0.785398163397448,0,0.785398163397448},',sprintf('\n\t'),...
    'xticklabels={{-$\pi$/2},{-$\pi$/4},{0},{$\pi$/4}},',sprintf('\n\t'),...
    'xlabel={$\prop_1$ (shoulder angle, rad)},',sprintf('\n\t'),...
    'ytick={0.785398163397448,1.5707963267949,2.35619449019234},',sprintf('\n\t'),...
    'yticklabels={{$\pi$/4},{$\pi$/2},{3$\pi$/4}},',sprintf('\n\t'),...
    'ylabel={$\prop_2$ (elbow angle, rad)},',sprintf('\n\t'),...
    'xmin=',num2str(S(1,1)),', xmax=',num2str(S(2,1)),', ',...
    'ymin=',num2str(S(1,2)),', ymax=',num2str(S(2,2)),', ',sprintf('\n\t'),...
    'width = ',num2str(S(2,1)-S(1,1)),'cm,',sprintf('\n\t'),...
  	'height = ',num2str(S(2,2)-S(1,2)),'cm,',sprintf('\n\t'),...
    'scale = ',num2str(scale),',',sprintf('\n'),...
    ']',sprintf('\n'),'\end{axis}',sprintf('\n\n')];



% loop through workspace
decoders = {'opt','EFH'};
for i = 1:Nx
    for j = 1:Ny
        
        % a single hand location, repeated Nexamples times
        thisTh = repmat([s1(i),s2(j)],[Nexamples,1]);
        thisX = FK2link(thisTh,params.roboparams,1);
        eStats = mastertest(wts,params,'numexamples',Nexamples,...
            'stimuli',cat(3,thisX,thisTh));
        
        % write the EFH and opt decoder error ellipses
        for iDecoder = 1:length(decoders)    
            ind = strcmp({eStats.tags(:).name},decoders{iDecoder});
            outtxt = [outtxt,tikzEllipse(...
                (eStats.Xpct(:,ind) + [s1(i),s2(j)]'),...
                chol(confScale*eStats.Cvrn(:,:,ind))',...
                [decoders{iDecoder},'clr'],0.5)];
        end
    end
end


% write file
titleStr = ['2DconditionalErrorStats-',params.NS,'-',date];
scaleStr = ['[x=',num2str(scale),'cm,y=',num2str(scale),'cm]'];
tikzWrite(outtxt,titleStr,scaleStr);


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