function states = sampler(means,FXN,params)
% SAMPLER   Sampler for Deep Belief Nets
%   SAMPLER produces a single sample (vector) drawn from the provided
%   MEANS, with the distribution specified by FXN: 'Gaussian', 'Bernoulli',
%   'Poisson', or 'Binomial'.
%
%   NB that to sample from a binomial distribution, SAMPLER converts the
%   MEANS (back) into sufficient statistics, first.
%
%   NB that params.t is *the number of Bernoulli neurons*, except for 'GP',
%   where it's the number of Gaussian neurons!!

%-------------------------------------------------------------------------%
% Revised: 03/??/14
%   -replaced poissrnd with ignpoi, which is about 20% faster in training
%   your standard RBM
% Revised: 05/21/13
%   -added 'GB' (Gaussian/Bernoulli)
% Revised: 11/30/10
%   -changed names of FXNs to rv type from inverse-link-fxn name
% Adapted: 9/28/10
%   -from updown.m or something0
%   by JGM
%-------------------------------------------------------------------------%

switch FXN
    case 'Gaussian' % w/fixed covariance          suff. stat. = mu
        states = means + randn(size(means),'like',means);
    case 'Bernoulli'                            % suff. stat. = pi = mu
        states = double(means > rand(size(means),'like',means));
        % states = round(means);
    case 'BernoulliDropout'                     % suff. stat. = pi = mu
        resp = means > rand(size(means),'like',means);
        resp(resp==1) = resp(resp==1).*(0.08>rand(sum(sum(resp)),1));
        states = double(resp);
%         states = double(...
%             means > rand(size(means),'like',means).*...
%             0.5 > rand(size(means),'like',means));
        % states = round(means);
    case 'Poisson'                              % suff. stat. = lambda = mu
        if strcmp(params.machine,'domestica')
            states = arrayfun(@ignpoiScalar,means);
        else
            states = ignpoi(means);
        end
    case 'Binomial'                             % suff. stat. = pi = mu/n
        states = binornd2(params.nexperiments,means/params.nexperiments);
    case 'PB'
        t = params.t;
        if strcmp(params.machine,'domestica')
            states = [arrayfun(@ignpoiScalar,means(:,1:end-t)),...
                double(means(:,end-t+1:end) > rand(size(means,1),t,'like',means))];
        else
            states = [ignpoi(means(:,1:end-t)),...
                double(means(:,end-t+1:end) > rand(size(means,1),t))];
        end
    case 'BP'
        t = params.t;
        if strcmp(params.machine,'domestica')
            states = [double(means(:,1:t) > rand(size(means,1),t,'like',means)),...
                arrayfun(@ignpoiScalar,means(:,(t+1):end))];
        else
            states = [double(means(:,1:t) > rand(size(means,1),t)),...
                ignpoi(means(:,(t+1):end))];
        end
    case 'GB'
        t = params.t;
        states = [means(:,1:end-t) + randn(size(means(:,1:end-t)),'like',means),...
            double(means(:,end-t+1:end) > rand(size(means,1),t,'like',means))];
    case 'BG'
        t = params.t;
        states = [double(means(:,1:t) > rand(size(means,1),t,'like',means)),...
            means(:,(t+1):end) + randn(size(means(:,(t+1):end)),'like',means)];
    case 'GP'
        % NB that here, t is the *number of Gaussians* rather than Bern.
        t = params.t;
        if strcmp(params.machine,'domestica')
            % states = [means(:,1:t) + 0.01*randn(size(means(:,1:t))),...
            states = [means(:,1:t) + randn(size(means(:,1:t)),'like',means),...
                arrayfun(@ignpoiScalar,means(:,(t+1):end))];
        else
            states = [means(:,1:t) + randn(size(means(:,1:t))),...
                ignpoi(means(:,(t+1):end))];
        end    
    case 'Dirac'
        states = means;
    otherwise
        error('unrecognized fxn');
end

                


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function r = binornd2(n,P)

r = sum(bsxfun(@gt,P,rand([size(P),n],'like',P)),3);
%%% r = sum(repmat(P,[1,1,n]) > rand([size(P),n]),3);

end
%-------------------------------------------------------------------------%
