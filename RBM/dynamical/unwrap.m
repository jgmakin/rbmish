function [x,steleport] = unwrap(x,posmin,posmax,steleport)
% UNWRAP takes src trajectories that are wrapped "like pac man" and unwraps them.
% It works by assuming that transitions are relatively close and you'll
% never see a jump as high as .7*range unless it is a wrap
%
% Inputs:
%   x: src, size ndim x T OR ncases x ndim x T and in world space.
%   posmin, posmax in world space
%
% Outputs:
%   x: unwrapped src
%   steleport: marks when unwrapping occured. Used to apply identical
%   unwrapping

%-------------------------------------------------------------------------%
% Revised: 06/7/13
%   -vectorized over ncases
% Revised: 06/6/13
%   -added steleport input option
% Revised: 06/3/13
%   -generalized to Ncases
% Created: 05/20/13
%   by BKD
%-------------------------------------------------------------------------%


% put in standard format: sz(x) = ncases x ndim x T
sz = size(x);
nsz = length(sz);
if nsz == 2
    x = reshape(x,[1,sz]);
end

[N,D,T] = size(x);

% if steleport exists, use the bounds it specifies
if exist('steleport','var')
    for d = 1:D
        posrange = posmax(d) - posmin(d);
        x(:,d,2:end) = x(:,d,2:end) - posrange.*cumsum(steleport(:,d,:),3);
    end
else
    % since steleport does not exist, attempt to find the jumps in the data
    steleport = NaN(N,D,T-1); %init
    
    for d = 1:D
        posrange = posmax(d) - posmin(d);
        pos = x(:,d,:);
        steleport(:,d,:) = (diff(pos,[],3) > .7*posrange) - (diff(pos,[],3) < -.7*posrange);
        x(:,d,2:end) = x(:,d,2:end) - posrange.*cumsum(steleport(:,d,:),3);
    end
end


end