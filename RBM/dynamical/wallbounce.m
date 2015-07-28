function zfuture = wallbounce(zfuture,zpast,smin,smax)
% Make trajectories bounce off walls.  The last two arguments are vectors
% containing the boundaries (a rectilinear workspace is assumed).

%-------------------------------------------------------------------------%
% Created: 05/11/13
%   by JGM
%-------------------------------------------------------------------------%


% check out all trajectories at once, or loop through them??

% params
m = length(smin);

% reflection matrices
R{1} = [-1 0; 0 1];
R{2} = [1 0; 0 -1];

% position, velocity
sfuture = zfuture(1:m);
spast = zpast(1:m);
vfuture = zfuture(m+1:2*m);

while any(sfuture < smin') || any(sfuture > smax')
    [scross,XingMod] = getCrossings(sfuture,spast,smin,smax);        
    sfuture = (sfuture - scross)*R{XingMod}' + scross;
    vfuture = vfuture*R{XingMod};
end

zfuture = [sfuture vfuture];

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [scross,crossingmod] = getCrossings(sfuture,spast,smin,smax)
%%% crossingmod can remain unassinged....

% init
d = sfuture - spast;
m = length(d);
inds = (1:m)';

% check for horizontal crossings, then for vertical crossings
for mod = 1:m
    
    if sfuture(inds(1)) > smax(inds(1))
        
        scross = getCross(smax,spast,d,inds);
        if (scross(inds(2))<=smax(inds(2)))&&(scross(inds(2))>=smin(inds(2)))
            crossingmod = inds(1);
        end
        
    elseif sfuture(inds(1)) < smin(inds(1))
        
        scross = getCross(smin,spast,d,inds);
        if (scross(inds(2))<=smax(inds(2)))&&(scross(inds(2))>=smin(inds(2)))
            crossingmod = inds(1);
        end
    end
    
    inds = circshift(inds,1);
    
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function scross = getCross(slim,spast,d,inds)
% this function will only work for m == 2

% xcross = (ycross - y1)*(x2-x1)/(y2-y1) + x1

scross(inds(1)) = slim(inds(1));
scross(inds(2)) = (slim(inds(1)) - spast(inds(1)))*...
    d(inds(2))/d(inds(1)) + spast(inds(2));

end
%-------------------------------------------------------------------------%