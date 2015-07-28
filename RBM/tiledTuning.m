function F = tiledTuning(X,g,inputUnitType,params)
% tiledTuning   A population of smoothly tiled tuning curves
%
% USAGE:
%   F = tiledTuning(X,g,params);
%
% Produces a matrix F of size (Nexamples x Nneurons) when given a matrix of
% inputs X (Nexamples x Ndims), a vector of gains g (Nexamples x 1), and
% the params structure.
%
% NB: main use is in PPCencode, but it can be used on its own for a
% "noiseless responses."

%-------------------------------------------------------------------------%
% Revised: 08/07/14
%   -replaced machine ('domestica') check with data class check
% Revised: 07/29/14
%   -added gpu stuff
% Revised: 07/24/14
%   -rewrote from scratch to use binary singleton expansion, tensorOp.m
%   -incorporated Ndims = 1, 2 to one case
%   -generalized to any Ndims
% Adapted: 06/21/14
%   -from GTrespfxn (version history below)
% Revised: 05/06/14
%   -vectorized inputs and outputs! (avoids parfor loop)
%   -renamed from respfxn to GTrespfxn
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% params
Info = inv(params.C);
N = params.N;
gridsize = params.gridsize;
[Nexamples, Ndims] = size(X);
Nunits = N^Ndims;

% build a lattice across Ndims with points at the PDs
PDs = linspace(0,gridsize,N)*ones(1,1,'like',X);
lattices1Dcell = mat2cell(PDs(ones(Ndims,1,'like',X),:),...
    ones(Ndims,1),N);
lattice = ndgrid(lattices1Dcell{:});
latticePDs = NaN(Nunits,1,Ndims,'like',X);
for iDim = 1:Ndims
    thisLattice = shiftdim(lattice,iDim-1);
    latticePDs(:,1,iDim) = thisLattice(:);
end
Infotensor = reshape(Info,[Ndims,1,Ndims])*ones(1,1,'like',X);
Y = bsxfun(@minus,reshape(X,[1,Nexamples,Ndims]),latticePDs);
YtrS = tensorOp(Y,Infotensor(:,ones(Nexamples,1,'like',X),:));
quadraticForm = sum(YtrS.*Y,3)';

%%% may ultimately want to make the tuning curves a parameter....
switch inputUnitType
    case {'Bernoulli','Binomial'}
        % subtract from logit(g), send through the logistic function
        logit = @(z)(log(z./(1-z)));
        logistic = @(z)(1./(1 + exp(-z)));
        F = logistic(bsxfun(@minus,logit(g),quadraticForm/2));
    otherwise % Gaussian tuning
        % exponentiate to get a Gaussian, scale by gains
        F = bsxfun(@times,g,exp(-quadraticForm/2));
end

end


