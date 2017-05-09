function [R,Q] = getDataCoherences(S,Q,dstrb)
% getDataCoherences
%
% USAGES:
%   R = getDataCoherences(S,Q)
%
% The counterpart to getStimuliCoherences.m, this function is used by the
% model 'ECcoherences'.

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from generateData.m
%   by JGM
%-------------------------------------------------------------------------%


Y = Q.coh(Q.inds);
if strcmp(dstrb,'GammaFixedScale')
    R = log(Y);
else % must be Erlang/Gamma
    R = [log(Y),Y];
end

end
