function h = areashade(x,y,f,c,thr)
% AREASHADE Shades areas between a curve and a fixed threshold
%
% AREASHADE(X,Y,F) shades the area where Y is larger than some fixed value
% F using the color yellow.
%
% AREASHADE(X,Y,F,C) uses color C instead of yellow.
%
% AREASHADE(X,Y,F,C,TH) where TH = 'h' uses F as a high threshold and
% shades the are where Y<F. Default for TH is 'l', painting Y>F.
%
% H = AREASHADE(X,Y,F) returns a handle to the created patch.
%
% Example:
% x = [0:.5:20]; y = sin(x);
% figure, plot(x,y,'o-'), hold on
% areashade(x,y,1/sqrt(2),'r')
% areashade(x,y,-1/sqrt(2),'b','h');
% plot(xlim,1/sqrt(2)*[1 1],'k')
% plot(xlim,-1/sqrt(2)*[1 1],'k')
%
% See also PATCH.

% aha, 24-mar-05, initial version

% default color yellow
if nargin<4 | isempty(c), c = 'y'; end
% define upper or lower threshold
if nargin<5 | isempty(thr), thr = 'l'; end

% flip to column vectors
if size(x,1) == 1
  x = x';
  y = y';
end

% add the crossover points (xf,f) to the data vectors to have a plotted
% line exactly on top of the shading
% from low to high
in = find(y(2:end)>f & y(1:end-1)<f);
m = (x(in+1)-x(in))./(y(in+1)-y(in));
xadd = x(in) + (f-y(in)).*m;
y = [y; f*ones(size(xadd))];
[x,si] = sort([x; xadd]);
y = y(si);
% from high to low
in = find(y(2:end)<f & y(1:end-1)>f);
m = (x(in+1)-x(in))./(y(in+1)-y(in));
xadd = x(in) + (f-y(in)).*m;
y = [y; f*ones(size(xadd))];
[x,si] = sort([x; xadd]);
y = y(si);




if thr=='h'
  % high threshold
  clip = y>=f;
elseif thr=='l'
  % low threshold
  clip = y<=f;
else
  error(['Supported threshold types: (h)igh or (l)ow: ' thr ' unknown!'])
end

y(clip) = f;
y = [f;y;f];
x = [x(1);x;x(end)];


h = patch(x,y,c,'edgecol','none');

% move patch into background
ax = get(h,'parent');
ch = get(ax,'Children');
try
  set(ax,'Children',[ ch(h~=ch) ch(h==ch)]);
catch
  set(ax,'Children',[ ch(h~=ch); ch(h==ch)]);
end

if nargout == 0
  clear
end

