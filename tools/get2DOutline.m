function s0 = get2DOutline(smin,smax,M)
%   s0 = get2DOutline(smin,smax,M)
%
%   get2DOutline returns a (fat) matrix (size 2xM) s0 that contains the 
%   outline of the square bounded by the 2-vectors smin and smax.

%-------------------------------------------------------------------------%
% Revised: 07/02/14
%   -transposed the output to be (N x 2)
%   -added "round"
% Cribbed: 01/09/14
%   from trajectoryGen.m
%   by JGM
%-------------------------------------------------------------------------%


srange = (smax - smin);
rat = srange(1)/srange(2);
MM = round(M*rat);
s0(:,1) = [linspace(smin(1),smax(1),MM)'; smax(1)*ones(M,1);...
    linspace(smax(1),smin(1),MM)'; smin(1)*ones(M,1)];
s0(:,2) = [smin(2)*ones(MM,1); linspace(smin(2),smax(2),M)';...
    smax(2)*ones(MM,1); linspace(smax(2),smin(2),M)'];


end