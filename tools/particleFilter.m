function [Xtu,W,XNtrpY] = particleFilter(Y,x0,transitionSampler,...
    emissionProb,USELOGS)
% particleFilter    For HMM-like models
%
% For examples of usage, see PFmaster.m in <code>/tools
%
% NB: if USELOGS = 1, the sampling weights W are log probabilities, not
% probabilities.

%-------------------------------------------------------------------------%
% Revised: 08/30/16
%   -changed output matrix X from "measurement update" to "time update"
%   -added weights as output argument
% Revised: 08/28/16
%   -moved "time update" (advancing sample along backbone) from the end to
%   the beginning of the loop, since otherwise the last X is at T+1 and it
%   is never measurement updated
% Created: 07/13/16
%   by JGM
%-------------------------------------------------------------------------%

% check arguments
if nargin < 5
    USELOGS = 0;
    if nargin < 4
        error('particleFilter.m needs more arguments -- jgm');
    end
end
T = size(Y,2);
[Nstates,Nparticles] = size(x0);
SLOWPROC = 1;

% open parallel pool
if Nparticles > 40000
    HADBEENCLOSED = isempty(gcp('nocreate'));
    if HADBEENCLOSED, pool = parpool(4); end
    catSmplFunc = @(w)(categorsmplPP(w,Nparticles));
    catLogSmplFunc = @(w)(categorlogsmplPP(w,Nparticles));
    %%% categorsmpl.m creates an enormous matrix....
else
    catSmplFunc = @(w)(categorsmpl(w,Nparticles,'IndexBased'));
    catLogSmplFunc = @(w)(categorlogsmpl(w,Nparticles));
end

% malloc/init
Xtu = NaN(Nstates,Nparticles,T,'like',Y);
W = NaN(Nparticles,T,'like',Y);
XNtrpY = 0;

if SLOWPROC, tic; fprintf('\n'); end
for t = 1:T

    % advance samples along the backbone, by sampling
    if t==1
        Xtu(:,:,t) = x0;
    else
        Xtu(:,:,t) = transitionSampler(Xmu,t-1);
    end
    
    % for numerical stability, you may want to use log-probabilities
    if USELOGS

        % construct the sampling weights from the likelihoods
        logpYgivenX = emissionProb(Y(:,t),Xtu(:,:,t),t)';    % Nsamples x 1
        logpYunnrmlzd = logprobs2logsumprobs(logpYgivenX);
        W(:,t) = logpYgivenX - logpYunnrmlzd;

        % accumulate cross entropy
        XNtrpY = XNtrpY - (logpYunnrmlzd - log(Nparticles));
        
        % re-draw (w/replacement) Nparticles sample trajectories
        trajInds = catLogSmplFunc(W(:,t));
        
    else
        % construct the sampling weights from the likelihoods
        pYgivenX = emissionProb(Y(:,t),Xtu(:,:,t),t)';
        pYunnrmlzd = sum(pYgivenX);
        W(:,t) = pYgivenX/pYunnrmlzd;

        % accumulate the cross entropy
        XNtrpY = XNtrpY - (log(pYunnrmlzd) - log(Nparticles));
        
        % re-draw (w/replacement) Nparticles sample trajectories
        trajInds = catSmplFunc(W(:,t));
    end
    Xmu = Xtu(:,trajInds,t);
    if SLOWPROC, fprintf('.'); end
end
if SLOWPROC, toc; fprintf('\n'); else fprintf('.'); end
if Nparticles > 40000, if HADBEENCLOSED, delete(pool); end; end
    
% time average the cross entropy (sensible)
XNtrpY = XNtrpY/T;      

end
