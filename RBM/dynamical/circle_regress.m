function [beta0,beta1] = circle_regress(x,y,mbounds)
%
%x is in arbitrary units, y is in radians, m0 is slope guess (default =
%randn

%example:
%{
a = 1 + randn*.5;
phi_0 = randn;

x = linspace(0,1,100);

ybounds = [0 10];
y = mod(a*x + phi_0 + randn(size(x))*.5,2*pi) + ybounds(1);

plot(x,y,'.')
xlabel('x')
ylabel('y')

[beta0,beta1] = circle_regress(x,y,[],ybounds,[-1 3]);

yhat = mod(beta1*x + beta0,2*pi);
hold all
plot(x,yhat,'.')
hold off
axis tight

%}

% check for varargin [where is m0 used??]
%%% if (~exist('m0','var') || isempty(m0)), m0 = randn; end
%%% if (~exist('wrap','var') || isempty(wrap)), wrap = [0 2*pi]; end

% a function whose derivative is proportional to circular variance
R = @(a)(1-sqrt((sum(cos(y - a*x)))^2 + (sum(sin(y - a*x)))^2)/length(y));

% find best-fit slope
if exist('mbounds','var') && ~isempty(mbounds)
    beta1 = fminbnd(R,mbounds(1),mbounds(2));
else
    beta1 = fminunc(R,randn);
end

% solve for y intercept
beta0 = atan2(sum(sin(y - beta1*x)),sum(cos(y - beta1*x)));


end





