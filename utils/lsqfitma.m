% lsqfitma.m					last update = 1 May 95
% 
% M-file to calculate a "2-way" least squares fit.
%       x,y data are passed as vectors from another m-file.
%       The line is fit by MINIMIZING the normal deviates.
%       The equation of the line is:     y = mx + b
%       This line is called the MAJOR AXIS.
%       Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
%            re-written from Kermack & Haldane (1950) Biometrika v37: 30-41;
%            after a derivation by Pearson (1901) Phil. Mag. v2(6): 559-572.
%
%       Data returned are as follows:
%
%	[m,b,r,sm,sb,xbar,ybar] = lsqfitma(x,y)
%
%               m    =    slope
%               b    =    y-intercept
%               r    =    correlation coefficient
%               sm   =    standard deviation of the slope
%               sb   =    standard deviation of the y-intercept
%               xbar =    mean of x values
%               ybar =    mean of y values
%
%       Note that (xbar,ybar) is the centroid
 
function [m,b,r,sm,sb,xbar,ybar]=lsqfitma(x,y)

% Determine the size of the vector
 
n = length(x);
 
% Calculate sums and other re-used expressions
 
Sx = sum(x);
Sy = sum(y);
xbar = Sx/n;
ybar = Sy/n;
u = x - xbar;
v = y - ybar;
 
Suv = sum(u .* v);
Su2 = sum(u .^2);
Sv2 = sum(v .^2);
 
sigx = sqrt(Su2/(n-1));
sigy = sqrt(Sv2/(n-1));
 
% Calculate m, b, r, sm, and sb
 
m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv);
b = ybar - m * xbar;
r = Suv / sqrt(Su2 * Sv2);
 
sm = (m/r) * sqrt((1 - r^2)/n);
sb1 = (sigy - sigx * m)^2;
sb2 = (2 * sigx * sigy) + ((xbar * m * (1 + r))/r^2);
sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);
