function mark(s,range,rmax,PANE,color,params)
%-------------------------------------------------------------------------%
% Revised: 01/24/11
%   -reduced to world2grid and the switch loop
% Revised: 01/21/11
%   -so as not to overlap so much with estError
% Cribbed: 01/21/11
%   from estimatorStats.m
%   by JGM
%-------------------------------------------------------------------------%

N = params.N;
Ndims = params.Ndims;
n = 100;                            % # pts. in lines; chosen arbitrarily

sGrid = world2grid(s,range,params);

switch Ndims
    case 1
        plot(sGrid*ones(1,n) + N*PANE,linspace(0,rmax,n),color);
    case 2
        plot(linspace(0,N,n) + N*PANE,sGrid(1)*ones(1,n),color);
        plot(sGrid(2)*ones(1,n) + N*PANE,linspace(0,N,n),color);
    otherwise
        fprintf('you never worked out crosses for Ndims!=2 -- jgm\n');
end