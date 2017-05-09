function T = sampleT(Mu,dstrbs,nums,params)
% sampleT   Sample sufficient-statistic generator for deep belief nets
%
% USAGE:
%
%   T = sampleT(Mu,'Poisson',params);  
%
% Given an Nexamples x (Nvars*Nstats) moment parameterization, Mu, of the
% factorial exponential-family distribution specified by dstrbs, first draw
% samples, then assemble the Nexamples x (Nvars*Nstats) sufficient stats.
%
% For most distributions, Nstats = 1, and the sample vector is itself the 
% minimal sufficient statistic ("linear suff. stats.") and is returned as
% T.  For those distributions with Nstats>1, the first Nvars columns of Mu
% and T correspond to the first statistic, and so on.  For example, for the
% Erlang and Gamma distributions, Nstats = 2; and the left and right halves
% of Mu and of T are assumed to hold the first and second first moments
% and sufficient statistics, respectively.  
% 
% NB: certain distributions require special care:
%   Erlang/Gamma
%   I don't know how to sample these using the moment parameterization only
%   and converting back and forth between shape and scale (on the one hand)
%   and moments (on the other) is expensive.  So "Mu" is really Th here,
%   which means these should not be substituted for states in EFH learning!
%
%   Binomial
%   Although Mu is indeed the first moment (the mean), i.e. p*n, n is
%   itself assumed to be known (and stored as params.Ntrials).


%-------------------------------------------------------------------------%
% Revised: 03/03/16
%   -added code to process faster (in one step) the case where both of the
%   dstrbs are the same
% Revised: 02/23/16
%   -replaced all "two-distribution" calls ('BP', etc.) with loop over cell
%   array of dstrbs and array of numsUnits
%   -added numsUnits as in input
% Revised: 02/12/16
%   -renamed from sampler.m to sampleT.m
% Revised: 01/21/16
%   -renamed some variables, and changed the comments accordingly
% Revised: 12/14/15
%   -added case for Erlang/Gamma
% Revised: 03/??/14
%   -replaced poissrnd with ignpoi, which is about 20% faster in training
%   your standard RBM
% Revised: 05/21/13
%   -added 'GB' (Gaussian/Bernoulli)
% Revised: 11/30/10
%   -changed names of FXNs to rv type from inverse-link-fxn name
% Adapted: 9/28/10
%   -from updown.m or something
%   by JGM
%-------------------------------------------------------------------------%

% if the number of unique groups is greater than 1, loop through them
if onlyOneUniqueStr(dstrbs)
    T = mmntParams2Samples2SuffStats(Mu,dstrbs{1},params);
else
    endinds = cumsum(nums);
    startinds = [1, endinds(1:end-1)+1];
    for iGrp = 1:length(dstrbs)
        T(:,startinds(iGrp):endinds(iGrp))=mmntParams2Samples2SuffStats(...
            Mu(:,startinds(iGrp):endinds(iGrp)),dstrbs{iGrp},params);
    end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function T = mmntParams2Samples2SuffStats(Mu,dstrb,params)

switch dstrb
    case 'Bernoulli'
        T = double(Mu > rand(size(Mu),'like',Mu));
    case 'BernoulliDropout'
        resp = Mu > rand(size(Mu),'like',Mu);
        resp(resp==1) = resp(resp==1).*(0.08>rand(gather(sum(sum(resp))),1));
        T = double(resp);
        %%% replace 0.08 with params.g??  Or something else???
        
        % T = double(...
        % 	 Mu > rand(size(Mu),'like',Mu).*...
        %    0.5 > rand(size(Mu),'like',Mu));
        % T = round(Mu);
    case 'StandardNormal'
        T = Mu + randn(size(Mu),'like',Mu);
    case 'Poisson'
        if isa(Mu,'gpuArray')
            T = arrayfun(@ignpoiScalar,Mu);
        else
            T = ignpoi(Mu);
        end
    case 'Binomial'
        T = binornd2(params.Ntrials,Mu/params.Ntrials);
    case {'Erlang','Gamma'}
        Th = max(Mu,0.0001);
        Y = gamrnd(Th(:,1:end/2),Th(:,(end/2+1):end));
        T = [log(Y), Y];
    case 'GammaFixedScale'
        Mu = max(Mu,0.0001);
        Y = gamrnd(gather(Mu),params.scaleparams);
		%%% this is sad
        T = log(Y);
		if isa(Mu,'gpuArray'), T = gpuArray(T); end
    case 'Categorical'
        T = categorsmpl(Mu',1,'OneHot')';
    case 'Dirac'
        T = Mu;
    otherwise
        error('unrecognized distribution');
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function r = binornd2(n,P)

r = sum(P > rand([size(P),n],'like',P),3);
%%% r = sum(repmat(P,[1,1,n]) > rand([size(P),n]),3);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function FLAG = onlyOneUniqueStr(strcell)

d = strcell{1};
FLAG = true;
i=1;
while FLAG&&(i<length(strcell))
    FLAG = strcmp(strcell{i+1},d);
    i=i+1;
end

end
%-------------------------------------------------------------------------%




% (2) attempt to compute (efficiently) *moment* parameters for Erlang dstrb
% Mu1 = thisMu(:,1:end/2);
% Mu2 = thisMu(:,(end/2+1):end);
% shapeparams = round(1./(2*(log(Mu2) - Mu1)));
% % you neglect terms at least as small as 1/(12*k^2)
% scaleparams = Mu2./shapeparams;
% shapeparams = max(shapeparams,0.0001);
% scaleparams = max(scaleparams,0.0001);
% Y = gamrnd(shapeparams,scaleparams);
