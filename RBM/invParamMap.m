function Mu = invParamMap(X,W,b,dstrbs,nums,params)
% invParamMap   Inverse parameter map on affine function of inputs
%  
% USAGE:
%
%   Mu = invParamMap(X,W,b,dstrbs,nums,params) 
%
% First an affine function is applied to the inputs, 
%
% 	X*W + b =: Eta.
%
% Then the result is passed through the "canonical" inverse parameter map,
%
%   Mu = ginv(Eta),
%
% i.e., the map that ensures that, for the (factorial) distribution Dstrbs,
% Mu corresponds to the *moment* parameterization and Eta to the *natural*
% parameters.  This guarantees that the natural parameters are linear fxns
% of the input.
%
% NB the Erlang and Gamma distributions are exceptions!!  Transforming from
% scale and shape params into the moments for Y and log(Y) is too expensive
% especially given that they just need to be transformed back in sampleT.m
% (it's unclear how to sample from these dstrbs using just the moment
% parameterization).
% 
%   Dstrbs              std params      inverse parameter map
%   ----------------------------------------------------------------
%   StandardNormal      means           Mu = Eta   
%   Bernoulli           means           Mu = 1/(1+exp{-Eta})
%   Binomial            means           Mu = 1/(1+exp{-Eta})
%   Poisson             means           Mu = exp{Eta}
%   Erlang/Gamma        shape, scale    "Mu" = [Eta_L + 1, -1/Eta_R]
%   Categorical         pi_{1:N-1}      Mu = exp{Eta}/(1 + sum(exp{Eta}))
%   ----------------------------------------------------------------
% 
% NB that the case Binomial assumes a known n (number of experiments) and
% that the "standard" parameter is taken to be *the mean*, i.e. p*n, and
% not just p.
%
% One can also request different functions for different "units"; see the
% case statements below.
%
% The *rows* of X are data vectors; each column of W is the set of weights 
% from all the input units to a single output unit; and the (row) vector b 
% contains the biases for each output unit.

%-------------------------------------------------------------------------%
% Revised: 03/03/16
%   -added code to process faster (in one step) the case where both of the
%   dstrbs are the same
% Revised: 02/23/16
%   -replaced all "two-distribution" calls ('BP', etc.) with loop over cell
%   array of dstrbs and array of numsUnits
%   -added numsUnits as an input
% Revised: 02/11/16
%   -renamed from feedforward.m to invParamMap.m
% Revised: 01/21/16
%   -renamed some variables, and changed the comments accordingly
%   -added case 'Categorical'
% Revised: 12/14/15
%   -added case 'Erlang'/'Gamma'
% Revised: 12/??/10
%   -added 'Binomial' as an option
% Revised: 09/27/10
%   -changed if/else to switch/case
% Revised: 06/02/10
%   -changed OUTFXN to a string
% Created: 05/25/10
%   by JGM
%-------------------------------------------------------------------------%

% "broadcast" (see note below)
Eta = X*W + b;

% if the number of unique groups is greater than 1, loop through them
if onlyOneUniqueStr(dstrbs)
    Mu = ntrlParams2mmntParams(Eta,dstrbs{1},params);
else
    endinds = cumsum(nums);
    startinds = [1, endinds(1:end-1)+1];
    for iGrp = 1:length(dstrbs)
        Mu(:,startinds(iGrp):endinds(iGrp)) = ntrlParams2mmntParams(...
            Eta(:,startinds(iGrp):endinds(iGrp)),dstrbs{iGrp},params);
           
    end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function Mu = ntrlParams2mmntParams(Eta,dstrb,params)

switch dstrb
    case {'Bernoulli','BernoulliDropout'}       % logistic
        Mu = 1./(1 + exp(-Eta));
    case 'StandardNormal'                       % linear
        Mu = Eta;
    case 'Poisson'                              % exponential
        Mu = exp(Eta);
    case 'Binomial'                             % scaled logistic
        Mu = 1./(1 + exp(-Eta))*params.Ntrials;
    case {'Erlang','Gamma'}
        Mu = [Eta(:,1:end/2)+1,  -1./Eta(:,(end/2+1):end)];
    case {'GammaFixedScale'}
        Mu = Eta+1;
    case 'Categorical'
        Mu = exp(Eta)./(1+sum(exp(Eta),2));
    case 'Dirac'
        Mu = Eta;
    otherwise
        error('unrecognized output function -- jgm\n');
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function FLAG = onlyOneUniqueStr(strcell)

d = strcell{1};
FLAG = true;
i = 1;
while FLAG&&(i<length(strcell))
   FLAG = strcmp(strcell{i+1},d);
   i=i+1;
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% NOTE. For future reference, it's helpful to know that the different ways
% of adding a vector to each column of a matrix don't appreciably change
% the amount of time it takes to compute X*W + b--perhaps b/c the matrix
% multiplication dominates so thoroughly.  Here are the approx. times:
%
% Eta = bsxfun(@plus,X*W,b);              % 121 s
% Eta = X*W + repmat(b,[size(x,1),1]);    % 122 s
% Eta = X*W + b(ones(size(x,1),1),:);     % 125 s
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% What's the right deterministic approximation for the Erlang distribution?
% That is, what should play the role of means in the Gaussian, Bernoulli,
% Poisson, etc. EFHs?  
% 
% We can re-phrase this as the reverse of a more common question:  Rather 
% than ask "What statistic is an unbiased estimator for mu?" we ask "For 
% what parameter is the sample mean an unbiased estimator?"  The answer is,
% for any distribution, the mean.  Thus, for Y ~ Erlang(shape,scale), the
% sample mean (\sum_i y_i)/N is an unbiased estimator of 
%
%   mu1 = shape*scale.
%
% Now, we don't know how Xi := log(Y) is distributed, but we do know that 
% the sample mean (\sum_i \log y_i)/N is an unbiased estimator of the mean 
% of Xi.  By the law of the unconscious statistician, we also know that the
% mean of Xi under p(xi) is the same as the mean of log(Y) under p(y).  And
% there is a formula for this expectation when p(y) is the Erlang dstrb:
%
%   E[log(Y)] = F(shape) + log(scale),
%
% where F is the *digamma* function.  For integer-valued shaped params k,
%
%   F(k) = H_{k-1} - gamma      = \sum_{n=1}^{k-1} 1/n   - gamma
%
% with H_n the nth harmonic number and gamma the Euler-Mascheroni constant.
% Thus, for log(Y), with Y ~ Erlang(shape,scale), the sample mean 
% (\sum_i \log y_i)/N is an unbiased estimator of
%
%   mu2 = \sum_{n=1}^{k-1} 1/n   - gamma + log(scale).
%
% Now we'll need to recover shape and scale from these (first) moments, so
% we can draw samples from gamrnd; i.e., we need the inverse function from
% mu1 and mu2 to the shape and scale parameters.
%
% Unfortunately, it's not clear how to invert the digamma function or the
% harmonic numbers.
%-------------------------------------------------------------------------%


% (1) Normexpon
%     case 'NORMEXPON'
%         th = exp(z);
%         th = th./repmat(sum(th,2),1,size(th,2));
%
% (2) attempt to compute (efficiently) *moment* parameters for Erlang dstrb
% shapeparams = Eta(:,1:end/2)+1;
% scaleparams = -1./Eta(:,(end/2+1):end);
% harmfunc = @(zz)(sum(1./(1:zz)));
% Mu1 = arrayfun(harmfunc,shapeparams-1) - double(vpa(eulergamma)) + log(scaleparams);
% Mu2 = shapeparams.*scaleparams;
% Mu = [Mu1, Mu2];
