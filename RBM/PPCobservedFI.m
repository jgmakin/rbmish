function FI = PPCobservedFI(Z,z,range,params)
% EMPFI     Empirical Fisher information
%   EMPFI computes the empirical or "observed" Fished information in a PPC
%   code Z, assuming Gaussian tuning curves and Poisson noise, around the
%   location z (which, generally speaking, will be the empirical mean).
%                                                   
%   Since Z is in grid coordinates and z is in real-world units, the range
%   of z must be provided to translate between them.
%-------------------------------------------------------------------------%
% Created: 11/30/10
%   by JGM
%-------------------------------------------------------------------------%


% init params
N = params.N;
g = params.g;
C = params.C;
respLength = params.respLength;
gridsize0 = params.gridsize;
r = range(:,2) - range(:,1);

% rescale (if necessary)
tunCov = (diag(r(:)./[respLength; respLength])*sqrtm(C))^2;
gridsize = scalefxn([gridsize0 gridsize0],[0 0],[respLength respLength],...
    [0 0],r);
S = inv(tunCov);

% compute firing rates f for this stimulus, z
b1 = linspace(0,gridsize(1),N)';
b2 = linspace(0,gridsize(2),N)';

d = (S(1,1)*repmat((z(1) - b1).^2,1,N) +...
    (S(2,1)+S(1,2))*(z(1) - b1)*(z(2) - b2)' +...
    S(2,2)*repmat((z(2) - b2').^2,N,1));
f = g*exp(-d/2);

% compute empirical fisher information
P = zeros(2,N);                                 % precompute grid coords.
for i = 1:N                                     %  in world units for speed
    P(:,i) = grid2world([i;i],range,params);
end
M = Z.^2 -  2*Z.*f + f.*f;                      % precompute for speed
FI = 0;
for i = 1:size(Z,1)                             % loop through neurons
    for j = 1:size(Z,2)
     
        zpref = [P(1,i); P(2,j)];               % preferred direction i,j
        q = S*(zpref - z);                      % dfbar/dz
        FI = FI + M(i,j)*(q*q')/4;
        
    end
end


end
