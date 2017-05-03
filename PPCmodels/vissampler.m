function v = vissampler(M,roboparams,varargin)
% VISSAMPLER    Sample uniformly from vis space
%   Vis space is shaped very strangely, so it's no easy task to sample even
%   uniformly from it.  VISSAMPLER addresses this issue by sampling first
%   from three squares, in proportion to their area; and then uniformly
%   from whichever square was chosen.  Together, these cover the space (see
%   wrkspc3.m for a visualization).  Samples which are outside vis space
%   are rejected, and the process restarted.

%-------------------------------------------------------------------------%
% Revised: 02/16/12
%   -changed to generate more than one sample at once; and (therefore)
%   internalized the outofrange checking (involves recursing)
% Created: 04/18/11
%   by JGM
%-------------------------------------------------------------------------%

% print warnings?
if nargin > 2; VERBOSE = varargin{1}; else VERBOSE = 1; end

posmin = roboparams.posmin;
posmax = roboparams.posmax;

% make three squares
%%%%%% these should really be derived from the params!!
xmin(:,1) = [-32,-21];
xmax(:,1) = [-22,24];

xmin(:,2) = [-22,-24];
xmax(:,2) = [-6,-8];

xmin(:,3) = [-22,18];
xmax(:,3) = [-8,32];

% get their areas
xarea = zeros(3,1);
for iSquare = 1:3
    xarea(iSquare) = prod(xmax(:,iSquare) - xmin(:,iSquare));
end

% pick a square (proportional to its area)
r = sum(xarea)*rand;
% iSquare = -1*(r < 450) + 1*(r > (xarea(1) + xarea(2))) + 2;
iSquare = (r > (xarea(1) + xarea(2))) - (r < 450) + 2;
% no, really, this is right

% pick within that square
x = scalefxn(rand(M,2),[0;0],[1;1],xmin(:,iSquare),xmax(:,iSquare))';

% rotate into the right place
phi = -pi/3;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
v = R*x;

% retain only the ones that are within the range
inds = zeros(size(v,2),1);
for i = 1:size(v,2)
    inds(i) = ~outofrange(v(:,i),posmin,posmax,'Hand-Position',0);
end
v = v(:,logical(inds));
P = M - size(v,2);

if P > 0
    % recurse to get the remainder
    if VERBOSE; fprintf('rejecting bad vis samples...'); end
    v1 = vissampler(P,roboparams,VERBOSE);
    v = [v v1];
else
    if VERBOSE; fprintf('\n'); end
end


end