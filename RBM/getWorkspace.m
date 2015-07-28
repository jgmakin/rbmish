function wrkspc = getWorkspace(M,params)
% gets the borders of the workspace, usually for plotting purposes

%-------------------------------------------------------------------------%
% Revised: 01/09/14
%   -pulled out get2Doutline into its own function
% Cribbed: 01/23/13
%   -from drawnoisyhills.m
%   by JGM
%-------------------------------------------------------------------------%


% trace the outline of prop space
th0 = get2DOutline(params.thmin,params.thmax,M);
x0 = FK2link(th0,params,1);

% store
wrkspc{strcmp(params.mods,'Hand-Position')} = x0;
wrkspc{strcmp(params.mods,'Joint-Angle')} = th0;


end