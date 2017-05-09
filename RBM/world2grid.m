function zGrid = world2grid(z,range,params)
% WORLD2GRID    converts stimuli from world to grid coordinates
%   WORLD2GRID performs the common transformation of the stimulus from its
%   true (world) representation to its coordinates in the grid (1:N).

%-------------------------------------------------------------------------%
% Created: 11/29/10
%   by JGM
%-------------------------------------------------------------------------%

% extract useful params
N = params.N;
Ndims = params.Ndims;
patchmin = params.margin*ones(1,Ndims);
patchmax = patchmin + params.respLength;

% convert x: world -> patch -> grid
zPatch = scalefxn(z,range(:,1),range(:,2),patchmin,patchmax);
zGrid = scalefxn(zPatch,zeros(Ndims,1),params.gridsize*ones(Ndims,1),ones(Ndims,1),...
    N*ones(Ndims,1));

end