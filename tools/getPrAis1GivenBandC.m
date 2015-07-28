function PrAis1givenBC = getPrAis1GivenBandC(PrAis1,b,c,Bedges,Cedges)
% PrAis1givenBC     Get Pr[A=1|B,C]
% 
% USAGE:
%   PrAis1givenBC = getPrAis1GivenBandC(PrAis1,b,c,Bedges,Cedges)
%
% Not infrequently---e.g., in dealing with neural networks---one has a
% vector of probabilities, PrAis1, of a Bernoulli random variables coming 
% up 1, for vectors of samples, b and c, from continuous random variables B
% and C (resp.).  Thus PrAis1 is at a discrete set of point, which is not
% very useful.  Instead, one would like to discretize B and C, and compute
%
%   Pr[A=1|B=b_d,C=c_d],
%
% i.e. the probability of heads given the discretized realizations of B and
% C.  This function returns that probability as PrAis1givenBC.  One also
% needs to supply the buckets (or rather their edges, Bedges and Cedges)
% into which B and C will be discretized.

%-------------------------------------------------------------------------%
% Created: 10/31/14
%   by JGM
%-------------------------------------------------------------------------%


% Ns
Nsamples = length(PrAis1);
[~,Bbins] = histc(b,Bedges);
[~,Cbins] = histc(c,Cedges);


% p(a,b)
Bincidence = sparse(1:Nsamples,Bbins,1,Nsamples,length(Bedges)-1);
Cincidence = sparse(1:Nsamples,Cbins,1,Nsamples,length(Cedges)-1);
pBC = full(Bincidence'*Cincidence);                 % omitting normalizer

% compute p(A=1,b,c) and p(A=1|b,c)
BincidenceTimesPofAis1 = sparse(1:Nsamples,Bbins,PrAis1,Nsamples,length(Bedges)-1);
PrAis1BC = full(BincidenceTimesPofAis1'*Cincidence); % omitting normalizer
PrAis1givenBC = PrAis1BC./(pBC + eps);
%%% replace eps trick with NaNs??


end