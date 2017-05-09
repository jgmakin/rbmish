function posteriorCumulants = gaussPosteriorization(cumulants)
% gaussPosteriorization     Get posterior from Gaussian likelihoods/priors
%
% USAGE:
%   posteriorCumulants = gaussPosteriorization(cumulants)
%
% Given a structure cumulants with fields Xpct (Nexample x Ndims x Nmods)
% and Info (Nexample x Ndims x Ndims x Nmods), which contain "cumulants"
% of "Gaussian" (log-quadratic) likelihoods or the cumulants of Gaussian 
% priors, returns a similar structure with the cumulants of the (Gaussian)
% posterior distribution, posteriorCumulants.  
%
% (See cumulantNeutralize.m for more information on the inputs.)
% 
% NB that the structures contain the *inverse* of the second cumulant, i.e.
% the information matrix.
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Revised: 08/07/14
%   -eliminated final loop through examples, replacing with arrayfun
% Revised: 07/04/14
%   -cleaned up, made it work; much of the old material has migrated to
%   cumulantNeutralize.m
% Revised: 07/03/14
%   -matricized computations
% Revised: 07/02/14
%   -incorporated integrateGaussianLikelihoods.m
% Created: 06/30/14
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Nexamples,Ndims,Nmods] = size(cumulants.Xpct);

% weighted combination
Xpct = permute(reshape(cumulants.Xpct,[Nexamples,Ndims*Nmods]),[2,3,1]);
Info = permute(reshape(cumulants.Info,[Nexamples,Ndims,Ndims*Nmods]),[2,3,1]);
posteriorCumulants.Info = sum(cumulants.Info,4);
infoINT = permute(posteriorCumulants.Info,[2 3 1]);
M = tensorOp(Info,Xpct);

% normalize by the integrated covariance--without a loop
if isa(cumulants.Xpct,'gpuArray')
    infoINT = gather(infoINT);
    M = gather(M);
end
I = repmat(reshape(1:Ndims*Nexamples,Ndims,1,Nexamples),[1,Ndims,1]);
J = repmat(reshape(1:Ndims*Nexamples,1,Ndims,Nexamples),[Ndims,1,1]);
INFOSPARSE = sparse(I(:),J(:),infoINT(:));
Xpct = reshape(INFOSPARSE\M(:),[Ndims,Nexamples])';

% now remove all the NaNs
badXpctInds = isnan(Xpct(:,1)); %%% just check dim=1; others are the same
if any(badXpctInds)
    
    % tell us how many
    fprintf('warning: found %d bad posterior means; ',sum(badXpctInds));
    fprintf('replacing with guesses -- jgm\n');
    
    % replace with...other, good posterior means!
    goodXpcts = Xpct(~badXpctInds,:);
    randInds = ceil(rand(sum(badXpctInds),1)*sum(~badXpctInds));
    Xpct(badXpctInds,:) = goodXpcts(randInds,:);
end
posteriorCumulants.Xpct = Xpct;

end
