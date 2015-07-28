function [R,gridShifts] = centerNoisyHills(R,N)
% centerNoisyHills      Just like it says
% 
% USAGE:
%   [R,gridShifts] = centerNoisyHills(R,params)
% 
% Given a matrix R (Nexamples x N^Ndims) of "firing rates," and the number,
% N, of neurons per dimension, returns R of the same size but in which the
% response volume for each example has been circularly shifted so as to put
% the neuron with the highest firing rate at the center of the volume.  The
% function also returns the amount of rotation along each dimension, for
% each example, in the matrix gridShifts (Nexamples x Ndims).
% 
% The code is highly vectorized and therefore very fast, totally
% impenetrable---but correct (use the commented out imagesc code to check).

%-------------------------------------------------------------------------%
% Created: 09/05/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
[M, Nlattice] = size(R);
Ndims = log(Nlattice)/log(N);
Nvec = [N*ones(1,Ndims,'like',R) ones(1,1,'like',R)];

% get the centering shifts
[~,theMaxInds] = max(R,[],2);
[maxSubscriptsCell{1:Ndims}] = ind2sub(Nvec,theMaxInds);
gridShifts = round((N*ones(M,Ndims,'like',R)...
    + 1)/2 - cat(2,maxSubscriptsCell{:}));

% get shifting indices for each row and column
inds = bsxfun(@minus,(0:N-1)',shiftdim(gridShifts,-1));
inds = mod(inds,N);


% now shift, looping only across dimsions
R = reshape(R,[M, Nvec]);                       % M x N x N x N...
for iDim = 1:Ndims
    R = permute(R,[2,1,3:(Ndims+1)]);           % swap way 1 and 2
    foo = inds(:,:,iDim);                       % row inds ->
    foo = bsxfun(@plus,foo,(1:N:N*M));          %  linear inds
    
    for jDim = 2:Ndims
        % make vec that skips across (jdim)D arrays
        vec = 0:(N^(jDim-1)*M):(N^jDim)*(M-1);
        
        % augment indices of jDim+1 way of foo
        foo = bsxfun(@plus,foo,shiftdim(vec,1-jDim));
    end
    R = shiftdim(R(foo),1);
end
R = reshape(R,[M,Nlattice]);
        

% VISUALIZE to check:
% for Ndims = 1:
%   plot(R(ceil(M*rand),:));
% for Ndims = 2:
%   imagesc(squeeze(R(ceil(M*rand),:,:)));
% for Ndims > 2:
% imagesc(squeeze(sum(R(ceil(M*rand),:,:,:),ceil(Ndims*rand + 1))))
 


end
