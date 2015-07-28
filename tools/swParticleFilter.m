%-------------------------------------------------------------------------%
%% init particle filter

% params
nsamples = 40000;

% malloc
emissprob = zeros(nsamples,1);
unnrmlklhd = zeros(T,1);
emisslogprob = zeros(nsamples,1);
LL = zeros(T,1);



%% run particle filter

% init samples (the same way you initialized the model)
musample = zeros(NY,nsamples,NZ);
for z = 1:NZ
    musample(:,:,z) = repmat(mu0(:,:,z),1,nsamples) + ...
        chol(CovM(:,:,z))*randn(NY,nsamples);
end


 
HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED, pool = parpool(4); end
tic
USELOGS = 0;
for t = 1:T
    
    % p(y_t|mu_t) = sum_{z_t}[p(y_t|z_t,mu_t)p(z_t)]
    if USELOGS
        y = Y(:,t)';
        parfor i = 1:nsamples
            b = mvnlpdf(y,squeeze(musample(:,i,:))',CovY) + log(pZ);
            emisslogprob(i) = logprobs2logsumprobs(b); 
        end
        unnrmlLL = logprobs2logsumprobs(emisslogprob);
        w = emisslogprob - unnrmlLL;
        LL(t) = unnrmlLL - log(nsamples);
        
        % push the samples forward one time step:
        % p(mu_{t+1}|y{0:t}) = <p(mu_{t+1}|musample_t)>_p(mu_t|y{0:t-1})
        whichsample = categorlogsmplPP(w,nsamples);
    else
        y = Y(:,t)';
        parfor i = 1:nsamples
            emissprob(i) = mvnpdf(y,squeeze(musample(:,i,:))',CovY)'*pZ;
        end
        unnrmlklhd(t) = sum(emissprob);
        w = emissprob/unnrmlklhd(t);
        
        % push the samples forward one time step:
        % p(mu_{t+1}|y{0:t}) = <p(mu_{t+1}|musample_t)>_p(mu_t|y{0:t-1})
        whichsample = categorsmplPP(w,nsamples);
    end
    
    for i = 1:NZ
        musample(:,:,z) = musample(:,whichsample,z) +...
            chol(CovM(:,:,z))*randn(NY,nsamples);
    end
    
end
toc
if HADBEENCLOSED, delete(pool); end



if ~USELOGS
    LL = log(unnrmlklhd/nsamples);
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
% precisely, in proportion as it was likely to have produced that
% emission.  The other mus will just evolve according to the set of m
% samples, i.e. without selecting a new multiset of samples---or more
% precisely, changing from the old samples to a new multiset of them only
% inasmuch as we believe the emission to have been generated from this mu.

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