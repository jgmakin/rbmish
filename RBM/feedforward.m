function y = feedforward(x,w,b,OUTFXN,params)
% FEEDFORWARD   Neural-network activation function
%   FEEDFORWARD(X,W,B) takes an affine function of data X with weights W,
%   offset by bias B, and passes the result through an output function:
%
%   linear:         Y = Z                           Gaussian
%   sigmoidal:      Y = 1/(1+exp{-Z})               Bernoulli/Binomial
%   exponential:    Y = exp(Z)                      Poisson
%
% One can also request different functions for different "units" (see case
% statements below).
%
% OUTFXN = 'Gaussian'/'Bernoulli'/'Poisson'/'Binomial' specifies this
% function by choosing the distribution's corresponding inverse canonical
% link.  The result---the sufficient statistic---is transformed into a mean
% if it isn't already.  THIS IS ONLY NECESSARY FOR 'Binomial'.
%
% The *rows* of X are data vectors; each column of W is the set of weights 
% from all the input (visible) units to a single output (hidden) unit; and
% the (row) vector B contains the biases for each output (hidden) unit.

%-------------------------------------------------------------------------%
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
z = bsxfun(@plus,x*w,b);

switch OUTFXN
    case {'Bernoulli','BernoulliDropout'}       % sigmoid
        y = 1./(1 + exp(-z));
    case 'Gaussian'                             % linear
        y = z;
    case 'Poisson'                              % exponential
        y = exp(z);
%     case 'NORMEXPON'
%         y = exp(z);
%         y = y./repmat(sum(y,2),1,size(y,2));
    case 'Binomial'
        y = 1./(1 + exp(-z))*params.nexperiments;
    case 'PB'
        t = params.t;
        y = [exp(z(:,1:end-t)), 1./(1 + exp(-z(:,end-t+1:end)))];
    case 'BP'
        t = params.t;
        y = [1./(1 + exp(-z(:,1:t))), exp(z(:,(t+1):end))];
    case 'GB'
        t = params.t;
        y = [z(:,1:end-t), 1./(1 + exp(-z(:,end-t+1:end)))];
    case 'GP'
        t = params.t;
        y = [z(:,1:t), exp(z(:,(t+1):end))];
    otherwise
        error('unrecognized output function -- jgm\n');
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% NOTE. For future reference, it's helpful to know that the different ways
% of adding a vector to each column of a matrix don't appreciably change
% the amount of time it takes to compute x*w + b--perhaps b/c the matrix
% multiplication dominates so thoroughly.  Here are the approx. times:
%
% z = bsxfun(@plus,x*w,b);              % 121 s
% z = x*w + repmat(b,[size(x,1),1]);    % 122 s
% z = x*w + b(ones(size(x,1),1),:);     % 125 s
%-------------------------------------------------------------------------%
