function z = decode(Z,range,params,MISC)
% decode
%   USAGE: 
%   z = decode(Z,range,params,MISC)
%
% NB: THIS FUNCTION SHOULD ALSO BE AVOIDE IN FUTURE: you have vectorized
% code for all of this

%-------------------------------------------------------------------------%
% Revised:
%   -put CoM and wrappeddecode in here, since they shouldn't be use
%   elsewhere: they've been superceded by tensorized versions.  
% Created: ??/??/?? (very early)
%   by JGM
%-------------------------------------------------------------------------%


% init
N = params.N;

% select decoding method
switch MISC
    case 'weighted'
        f = MISC{2};
        Zw = zeros(N,N);
        % weighted (by firing fxn) sum
        for i = 1:N
            for j = 1:N
                % input
                fmat = f(1+i-1:N+i-1,1+j-1:N+j-1);
                Zw(i,j) = sum(sum(fmat.*Z));
            end
        end
        [Z_a,Z_b] = find(Zw == max(max(Zw)));
        zhat = [Z_a Z_b];
    
    case 'CoM'                                  % alt.: center of mass
        zhat = CoM(Z);
        
    case 'torusCoM'
        zhat = wrappeddecode(Z,params);         % alt.: CoM on a torus
    
    case 'template'                             % alt.: template matching
        f = MISC{2};
        METRIC = 'Euclidean';
        [Z_a,Z_b] = templatematcher(Z,f,METRIC);
        zhat = [Z_a Z_b];
    otherwise
        error('unrecognized estimator -- jgm');
end

% rescale: grid -> patch -> world [see note]
z = grid2world(zhat,range,params);

end
% NB that you rescale into the whole grid, [0,gridsize], rather than the
% centered patch, [margin,resplength+margin].  That's correct, since you
% want it to be *theoretically* possible that CoM put the center of mass
% outside the range of sample thetas and positions.  Generally speaking, of
% course, that won't happen.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function cntr = CoM(A)
% CoM   Matrix center of mass
%   USAGE: cntr = CoM(A)
%
%   Compute the center of mass of a matrix A, where the weight of each
%   particle is given by an entry in A and its position is given by the
%   indices of that entry.
%
%   NB: The order of the cntr coordinates matches matlab's dimension
%   ordering (row, column,...)
%
%   NB: If the array order is 1 (i.e., if params.m, the dimension of the
%   stimulus, is 1), then CoM expects a *column* vector (cf. displayshape).
%-------------------------------------------------------------------------%
% Revised: 12/06/10 (happy b'day)
%   -rewrote to work for any size array, not just 2D (thereby making this
%       fxn even more opaque)
% Created: 07/5/10
%   by JGM
%-------------------------------------------------------------------------%

denom = sum(A(:));
arrayorder = sum(size(A)~=1);
cntr = zeros(arrayorder,1);

for i = 1:arrayorder
    if denom == 0
        cntr(i) = randsample(size(A,i),1);
        fprintf('warning: CoM detects a zero vector -- jgm\n');
    else
        cntr(i) = sum([1:size(A,i)]*shiftdim(A,i-1))/denom;
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function z = wrappeddecode(Z,params)
% For decoding on a torus!

%-------------------------------------------------------------------------%
% Retired: 07/07/14
%   -still in circulation, but this has been superceded by a function which
%   does most of the calculations matrix/tensor-wise, and lives in
%   GTPNsuffstats.m.
% Created: ??/??/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
N = params.N;
Ndims = params.Ndims;

% center around the maximal point
[~,theMaxInd] = max(Z(:));
[maxSubscriptsCell{1:Ndims}] = ind2sub(size(Z),theMaxInd);
maxSubscripts = [maxSubscriptsCell{:}];
gridShifts = round((N*ones(1,Ndims) + 1)/2 - maxSubscripts);
Zcentered = circshift(Z,gridShifts);

% get the center of mass of the *centered* response
zCntrdhat = CoM(Zcentered);

% shift back into place
z = zCntrdhat - gridShifts';


end
%-------------------------------------------------------------------------%