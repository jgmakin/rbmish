function [PrAis1givenBC,MI_AwithBC] =...
    getPrAis1givenBandC(PrAis1givenBCcont,b,c,Bedges,Cedges)
% getPrAis1givenBandC     Get Pr[A=1|B,C] and I[A;B,C]
% 
% USAGE:
%   PrAis1givenBC = getPrAis1GivenBandC(PrAis1GivenBCcont,b,c,Bedges,Cedges)
%
% Not infrequently---e.g., in dealing with neural networks---one has a
% vector of probabilities, PrAis1givenBCcont, of Bernoulli random variable
% coming up 1, for vectors of samples, b and c, from continuous random 
% variables B and C (resp.).  Thus PrAis1givenBCcont is at a discrete set 
% of point, which is not very useful.  Instead, one would like to 
% discretize B and C, and compute
%
%   Pr[A=1|B=b_d,C=c_d],
%
% i.e. the probability of heads given the *discretized* realizations of B 
% and C.  This function returns that probability as PrAis1givenBC. One also
% needs to supply the buckets (or rather their edges, Bedges and Cedges)
% into which B and C will be discretized.
% 
% The second output argument, MI_AwithBC, is the mutual information:
%
%   I(A; B,C).
%
% NB: PrAis1givenBCcont must have size Nsamples x Nunits (Nunits can be 1).  
% The output PrAis1givenBC has size NBbins x NCbins x Nunits.  Note that 
% NBbins (e.g.) is length(NCedges) - 1.

%-------------------------------------------------------------------------%
% Revised: 07/29/16
%   -added mutual-information output
% Created: 10/31/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
[Nsamples,Nunits] = size(PrAis1givenBCcont);
if length(b) ~= Nsamples
    error('Nsamples in b and PrAis1 don''t match -- jgm');
end
[~,Bbins] = histc(b,Bedges);
[~,Cbins] = histc(c,Cedges);

% Pr(A,B)
Bincidence = sparse(1:Nsamples,Bbins,1,Nsamples,length(Bedges)-1);
Cincidence = sparse(1:Nsamples,Cbins,1,Nsamples,length(Cedges)-1);
PrBC = full(Bincidence'*Cincidence);                 % omitting normalizer

% malloc
PrAis1givenBC = zeros([size(PrBC),Nunits],'like',PrAis1givenBCcont);

% Pr(A=1,B,C), Pr(A=1|B,C)
for iUnit = 1:Nunits
    BincidenceTimesPofAis1 = sparse(1:Nsamples,Bbins,PrAis1givenBCcont(:,iUnit),...
        Nsamples,length(Bedges)-1);
    PrAis1BC = full(BincidenceTimesPofAis1'*Cincidence); % omitting normalizer
    PrAis1givenBC(:,:,iUnit) = PrAis1BC./(PrBC + eps);
    %%% replace eps trick with NaNs??
end



% mutual information between A (on the one hand) and B and C (on the other)
if nargout > 1
    
    % P(A=0|B,C), P(A=1|B,C)
    PrAis0givenBC = 1 - PrAis1givenBC;
    PrAis1givenBC = PrAis1givenBC+eps;
   
    % H[A|B,C]
    HAis0givenBC = sum(sum((-log2(PrAis0givenBC).*PrAis0givenBC).*PrBC,1),2);
    HAis1givenBC = sum(sum((-log2(PrAis1givenBC).*PrAis1givenBC).*PrBC,1),2);
    HAgivenBC = (squeeze(HAis0givenBC)' + squeeze(HAis1givenBC)')/Nsamples;
   
    % H[A]
    PrA = [1-sum(PrAis1givenBCcont,1)/Nsamples; sum(PrAis1givenBCcont,1)/Nsamples];
    HA = sum(-log2(PrA).*PrA,1);
    
    % I[A;B,C] = H[A] - H[A|B,C]
    MI_AwithBC = HA - HAgivenBC;
end


end
