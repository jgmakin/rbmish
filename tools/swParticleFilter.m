%-------------------------------------------------------------------------%
% Replaced: 07/14/16
%   -by spikesorted.m using particleFilter.m
% Revised: 07/13/16
%   -replaced loops over samples with bsxfun and loops over NZ (and even
%   these can probably be eliminated)
% Created: ??/??/12
%   by JGM
%-------------------------------------------------------------------------%

%% init particle filter

% params
Nsmpls = 40000;

% malloc
condEmissProbs = zeros(Nsmpls,NZ);
condEmissLogProbs = zeros(Nsmpls,NZ);
unnrmLklhd = zeros(T,1);
LL = zeros(T,1);



%% run particle filter

% init samples (the same way you initialized the model)
Msmpls = zeros(NY,Nsmpls,NZ);
for z = 1:NZ
    Msmpls(:,:,z) = mu0(:,:,z) + chol(CovM(:,:,z))'*randn(NY,Nsmpls);
end

 
HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED, pool = parpool(4); end
tic
USELOGS = 0; % more accurate (presumably), but takes about 5.5x time
for t = 1:T
    
    % observations at time t
    y = Y(:,t)';
    
    
    if USELOGS
        % compute log likelihoods for each discrete state (for all samples)
        for z = 1:NZ
            condEmissLogProbs(:,z) = mvnlpdf(y,Msmpls(:,:,z)',CovY(:,:,z));
        end
        
        % marginalize out Z
        % p(y_t|mu_t) = sum_{z_t}[p(y_t|z_t,mu_t)p(z_t)]
        emissLogProbs = logprobs2logsumprobs(condEmissLogProbs' + log(pZ))';
        
        % construct sampling weights
        unnrmlLL = logprobs2logsumprobs(emissLogProbs);
        w = emissLogProbs - unnrmlLL;
        
        % advance along the backbone: discrete-state samples
        % p(z_{t+1}|y{0:t}) = <p(z_{t+1}|z_t)>_p(z_t|y{0:t})
        zSmpls = categorlogsmplPP(w,Nsmpls);
        %%% can this be sped up?
        
        % store the log-likelihood of the data under the model
        LL(t) = unnrmlLL - log(Nsmpls);
    else
        % compute likelihoods for each discrete state (for all samples)
        for z = 1:NZ
            condEmissProbs(:,z) = mvnpdf(y,Msmpls(:,:,z)',CovY(:,:,z));
        end
        
        % marginalize out Z
        % p(y_t|mu_t) = sum_{z_t}[p(y_t|z_t,mu_t)p(z_t)]
        emissProbs = condEmissProbs*pZ;
        
        % construct sampling weights
        unnrmLklhd(t) = sum(emissProbs);
        w = emissProbs/unnrmLklhd(t);
        
        % advance along the backbone: discrete-state samples
        % p(z_{t+1}|y{0:t}) = <p(z_{t+1}|z_t)>_p(z_t|y{0:t})
        zSmpls = categorsmplPP(w,Nsmpls);
        %%% categorsmpl.m creates an enormous matrix....
    end
    
    % advance along the backbone: cluster-mean samples
    % p(mu_{t+1}|y{0:t}) = <p(mu_{t+1}|mu_t)>_p(mu_t|y{0:t})
    for z = 1:NZ
        Msmpls(:,:,z) = Msmpls(:,zSmpls,z) + chol(CovM(:,:,z))'*randn(NY,Nsmpls);
    end
    
end
toc
if HADBEENCLOSED, delete(pool); end


% compute the log-likelihood of the data under the model
if ~USELOGS
    LL = log(unnrmLklhd/Nsmpls);
end
%-------------------------------------------------------------------------%

% may need a lot of samples b/c of C.o.D.: you need to sample from the
% product of:
%   p(a_{t+1}|a_t)p(b_{t+1}|b_t)p(c_{t+1}|c_t)p(d_{t+1}|d_t)...,
% where (a,b,c,d...) = the set of mu's.  These are based on samples of a_t,
% b_t, c_t, d_t,....  With four letters, you've got 4 NY-dimensional
% spaces, which is pretty big....



% the idea will almost certainly be that you just let the "unobserved"
% data evolve without updating---i.e., you'll fix m samples *per* mu; then
% you'll update to p(mu_{t+1}|y{0:y}), using the transition probabilities,
% only for the mu which was likely to have produced the emission---or more
% precisely, in proportion as it was likely to have produced that emission.
% The other mus will just evolve according to the set of m samples, i.e.,
% without selecting a new multiset of samples---or more precisely, changing
% from the old samples to a new multiset of them only inasmuch as we 
% believe the emission to have been generated from this mu.

% The intuition here is that if e.g. you've got some good samples of mu1
% (call it a) in the same vector as some bad samples of mu2 (call it b),
% then you want to be able to exploit that goodness w/o suffering from the
% badness: you want to be able to propagate along the a's with high
% probability, while throwing out the bad b's that have been "paired with
% it."



%% exact log-likelihoods

Lcum = 0;
LLreal = zeros(T,1);
for i = 1:T
    L = mvnlpdf(Y(:,i)',squeeze(mu0)',CovY+i*CovM) + log(pZ);
    Lcum = repmat(Lcum(:),1,NZ) + repmat(L',numel(Lcum),1);
    LLreal(i) = logprobs2logsumprobs(Lcum(:));
end    
logprobs2logsumprobs(Lcum(:));

% pf estimate
cumsum(LL)
LLreal

    
    
    
    
% to do:
% (1) check what's costing you with the profiler
% (2) plot out results w/different numbers of particles