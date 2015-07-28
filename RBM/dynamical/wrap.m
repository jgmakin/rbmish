function out = wrap(in,posmin, posmax)
% WRAP takes an unwrapped src and makes it wrapped (like pacman)
% reverses UNWRAP

% inputs:
%   in: usrc, size DxT and in world space
%   posmin, posmax: 1xD and in world space

% outputs:
%   out: src

D = size(in,1);

for d = 1:D
    
    range = posmax(d) - posmin(d);
    out = mod(in(d,:)-posmin(d),range) + posmin(d);
    
end