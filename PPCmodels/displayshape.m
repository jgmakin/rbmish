function topogrids = displayshape(d,params)
% DISPLAYSHAPE  Topographically organize
%   USAGE:
%       topogrids = displayshape(d,params);
%
%   DISPLAYSHAPE turns vectors d into params.Nmods multidimensional arrays with
%   params.m dimensions each.

%-------------------------------------------------------------------------%
% Revised: 12/17/10
%   -rewrote completely to accomodate #modalities != 2
% Revised: 12/08/10
%   -removed non-square-grid stuff and added accomodations for non-2D
%   arrays.
% Created: 06/22/10
%   by JGM
%-------------------------------------------------------------------------%

Nmods = length(params.mods);
Ndims = params.Ndims;
N = params.N;

% put each of Nmods modalities in its own column
D = reshape(d,length(d)/Nmods,Nmods);

% reshape according to the dimension Ndims of the stimulus
topogrids = cell(1,Nmods);
for i = 1:Nmods
    topogrids{i} = reshape(D(:,i),[N*ones(1,Ndims),1]);
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [a,b] = sqrtapprox(P)

a = floor(sqrt(P));
while P/a - floor(P/a) > 0
    a = a - 1;  
end
b = P/a;

end
%-------------------------------------------------------------------------%

% function f = FermatFactor(N)                % N should be odd
% a = ceil(sqrt(N));
% b2 = a*a - N;
% while (sqrt(b2) - floor(sqrt(b2))) > 0
%     a = a + 1;
%     b2 = a*a - N;                           % equiv.: b2 ? b2 + 2*a - 1
% end
% f = a - sqrt(b2);                           % or a + sqrt(b2)
% end