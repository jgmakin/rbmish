function shatL = decoder(d,params)
% DECODER   Decode a vector of input populations
% DECODER decodes an entire vector of observed activities (all populations)
% NB that this runs only on single data vectors; run in a loop to get all
% the decodings for a whole data set D0.
%
% NB: This function will also remove the mark that DATAGENPP puts into
% decoupled data.

%-------------------------------------------------------------------------%
% Created: 02/15/12
%   by JGM
%-------------------------------------------------------------------------%


% init
Nmods = params.Nmods;
Ndims = params.Ndims;
smin = params.smin;
smax = params.smax;
shatL = zeros(Ndims,Nmods);

% decode
%%% hacky
% d(params.N) = 0;
%%%
T = displayshape(d,params);
for i = 1:Nmods
    shatL(:,i) = decode(T{i},[smin(:,i) smax(:,i)],params,'CoM');
end

end
