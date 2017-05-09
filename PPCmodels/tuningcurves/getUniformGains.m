function G = getUniformGains(M,params)
% generates a uniform MxM grid of gains

%-------------------------------------------------------------------------%
% Created: 07/27/12
%   by JGM
%-------------------------------------------------------------------------%

%%% params.g is deprecated....
foo = linspace(0,params.g*2,M);
[gains1 gains2] = meshgrid(foo,foo);
G = [gains1(:),gains2(:)];
%%% is this really the right range---gains in [0,2*params.g] ???

end