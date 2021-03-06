function samples = categorsmpl(pmf,Ntrials,representationType)
% categorsmpl   Sampler from a categorical distribution
%   categorsmpl(pmf,Ntrials,representationType) draws Ntrials samples from 
%   a categorical distribution with probabilities pmf, a matrix of size 
%   Ncats x Ndice.
%
% NB that pmf must sum to one across *rows*.
%
% The output form is determined by representationType, which can be either
% 'OneHot', or something else (suggested string is 'IndexBased')--in which 
% case the index of the category is returned.
% 
% NB the shape of the returned matrix of samples: For 'OneHot', it has size
%
%   Ncats x Ndice x Ntrials
%
% whereas for non-OneHot representations it has size
%
%   Ntrials x Ndice.

%-------------------------------------------------------------------------%
% Revised: 07/12/16
%   -re-wrote (entirely) w/bsxfun instead of repmat, and allowing matrix
%   inputs for multiple trials
%   -added flag representationType to toggle one-hot or index-based.
% Created: 03/16/12
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Ncats,Ndice] = size(pmf);

% don't create too large of a tensor 
if Ncats*Ndice*Ntrials > 100000000
    fprintf('Ncats*Ndice*Ntrials > 100,000,000; recursing...   ');
    fprintf('(OneHot representation not advised)\n');
    
    % break into two pieces
    samplesA = categorsmpl(pmf,floor(Ntrials/2),representationType);
    samplesB = categorsmpl(pmf,Ntrials - floor(Ntrials/2),representationType);
    if ~strcmp(representationType,'OneHot')
        samples = cat(1,samplesA,samplesB);
    else
        % but this is probably too big
        samples = cat(3,samplesA,samplesB);
    end
    
else
    % the sampling procedure
    samples = diff([zeros(1,Ndice); cumsum(pmf)] > rand(1,Ndice,Ntrials));
    
    if ~strcmp(representationType,'OneHot')
        [samples, ~] =  find(samples);
        samples = reshape(samples,[Ndice,Ntrials])';
    end
    
end

end





