function [PrAis1givenB,MI_AwithB] = getPrAis1givenB(PrAis1givenBcont,b,Bedges)
% getPrAis1givenBandC     Get Pr[A=1|B] and I[A;B]
%
% USAGE:
%   PrAis1givenB = getPrAis1GivenB(PrAis1GivenBcont,b,Bedges)
%
% Not infrequently---e.g., in dealing with neural networks---one has a
% vector of probabilities, PrAis1givenBcont, of Bernoulli random variable
% coming up 1, for vectors of samples, b, from a continuous random variable
% B.  Thus PrAis1givenBcont is at a discrete set of points, which is not 
% very useful.  Instead, one would like to discretize B, and compute
%
%   Pr[A=1|B=b_d],
%
% i.e. the probability of heads given the *discretized* realizations of B. 
% This function returns that probability as PrAis1givenB. One also needs to
% supply the buckets (or rather their edges, Bedges) into which B will be 
% discretized.
%
% NB: getPrAis1givenB must have size Nsamples x Nunits! (Nunits can be 1).  
% The output PrAis1givenBC has size NBbins x Nunits.  Note that NBbins is 
% length(NCedges) - 1.


%-------------------------------------------------------------------------%
% Revised: 07/29/16
%   -re-wrote to accommodate *matrix* inputs PrAis1givenBcont
% Cribbed: 07/29/16
%   from getOptimalLag.m
%   by JGM
%-------------------------------------------------------------------------%


% Ns
Nsamples = size(PrAis1givenBcont,1);
NBbins = length(Bedges)-1;

% Pr(B)
[Nb,Bbins] = histc(b,Bedges);
PrB = Nb(1:NBbins)/Nsamples + eps;  %%% just in case

% Pr(A,B)
Bincidence = sparse(1:Nsamples,Bbins,1,Nsamples,NBbins);
PrAis1B = (Bincidence'*PrAis1givenBcont)/Nsamples;
PrAis1B(PrAis1B==0) = PrAis1B(PrAis1B==0) + eps;

% Pr(A=1|B)
PrAis1givenB = PrAis1B./PrB;


% mutual information between A and B
if nargout > 1
    
    % Pr(A=0|B)
    PrAis0B = (Bincidence'*(1-PrAis1givenBcont))/Nsamples;
    PrAis0B(PrAis0B==0) = PrAis0B(PrAis0B==0) + eps;
   
    % Pr(A)
    PrAis0 = 1-sum(PrAis1givenBcont)/Nsamples; 
    PrAis1 = sum(PrAis1givenBcont)/Nsamples;
    
    % I(A;B), mutual information
    MI_AwithB = ...
        sum(PrAis0B.*(log2(PrAis0B) - log2(PrB*PrAis0))) +...
        sum(PrAis1B.*(log2(PrAis1B) - log2(PrB*PrAis1)));
end


end
