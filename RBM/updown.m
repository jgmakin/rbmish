function [err,Rout] = updown(Rin,wts,params,varargin)
% UPDOWN   Propagation through a DBN
%   UPDOWN(D0,WTS,PARAMS) propagates D0 through the DBN specified by WTS
%   and PARAMS, either (depending on your pov) all the way through the
%   "unfolded" network," or up and down the original one (tho' using
%   different biases for each direction).
%
%   You can additionally specify the propagation method ('samples','MLE',
%   'modes','means'), where the default is something strange (probably a
%   combination of different propagation methods for different units).
%
%   Its original use in fine-tuning the deep autoencoder, but it turns out
%   to be useful in a lot of places.
%
% NB: This function expects longdata!!!

%-------------------------------------------------------------------------%
% Revised: 12/12/13
%   -changed to use inputs Rin that are in longdata format only
% Revised: 05/21/13
%   -added options for 'GB'
% Revised: 04/26/12
%   -added options for 'BP'
% Revised: 09/14/11
%   -renamed updown.m (the old updown was changed into updownfast.m)
%   -added some options for mixed propagation methods
% Revised: 09/08/11
%   -swapped the ordering of the switch statements, etc. etc.
%   -added options for 'PB'
% Revised: 08/03/10
%   -add params as an argument and fixed the way FXN work---now you have to
%   make sure it's set up properly in setParams
% Revised: 5/26/10
%   -added dataout as an output
% Adapted: 5/25/10
%   from the Salakhutdinov/Hinton code
%   by JGM
%-------------------------------------------------------------------------%

% init
DISP = 0;
VERBOSE = 1;
numlayers = length(wts)/2+1; % length(params.numsUnits);
TYPE = 'default';

for i=1:length(varargin)
    switch varargin{i}
        case {'samples','MLE','modes','means','Nsamples','poissonmodes'}
            TYPE = varargin{i};
        case 'quiet'
            VERBOSE = 0;
    end
end
if VERBOSE
    fprintf('using method %s\n',TYPE);
end


% cycle through layers (twice as many as the number of RBMs)
states = Rin;
for layer = 1:length(wts)
    FXN = params.typeUnits{numlayers-abs(numlayers-layer-1)};
    % (really)
    
    means = feedforward(states,wts{layer}(1:end-1,:),...
        wts{layer}(end,:),FXN,params);
    if ~sum(strcmp(FXN,{'Gaussian','Bernoulli','Poisson','Binomial',...
            'PB','BP','GB','GP','BernoulliDropout'}))
        error('unrecognized NN fxn - jgm\n');
    end
    
    switch TYPE
        case 'samples'
            states = sampler(means,FXN,params);
        case 'means'
            switch FXN
                case 'BernoulliDropout'
                    states = feedforward(states,100*wts{layer}(1:end-1,:),...
                        wts{layer}(end,:),FXN,params);
                otherwise
                    states = means;
            end
        case 'modes'
            switch FXN
                case 'Gaussian'
                    states = means;
                case 'Poisson'
                    states = floor(means);
                case 'Bernoulli'
                    states = round(means);
                case 'Binomial'
                    n = params.nexperiments;
                    states = floor(means*(n+1)/n);
                case 'PB'
                    t = params.t;
                    states(:,1:end-t) = floor(means(:,1:end-t));
                    states(:,end-t+1:end) = double(means(:,end-t+1:end) > 1/2);
                case 'BP'
                    t = params.t;
                    states(:,1:t) = double(means(:,1:t) > 1/2);
                    states(:,(t+1):end) = floor(means(:,(t+1):end));
                case 'GB'
                    t = params.t;
                    states(:,1:end-t) = means(:,1:end-t);
                    states(:,end-t+1:end) = double(means(:,end-t+1:end) > 1/2);
                case 'GP'
                    t = params.t;
                    states(:,1:t) = means(:,1:t);
                    states(:,(t+1):end) = floor(means(:,(t+1):end));
            end
        case 'poissonmodes'                 % for "output decoding"
            switch FXN
                case 'Poisson'
                    states = floor(means);
                otherwise
                    states = means;
            end
        case 'Nsamples'                     % for "hidden layer decod."
            switch FXN
                case 'Bernoulli'
                    NNNN = params.smpls;
                    states = sampler(means,FXN,params);      % samples!
                    for i=1:NNNN-1
                        states = states + sampler(means,FXN,params);
                    end
                    
                    if VERBOSE
                        fprintf('using %i Bernoulli sample(s)...\n',NNNN);
                    end
                    states = states/NNNN;
                    
                otherwise
                    states = means;
            end
        otherwise
            switch FXN
                case 'Poisson'
                    states = round(means); % means - 0.5;
                    if VERBOSE
                        fprintf('*rounding* Poisson means...\n');
                    end
                case 'PB'
                    t = params.t;
                    sts1 = floor(means(:,1:end-t));
                    fprintf('using modes for Poisson neurons...\n');
                    sts2 = means(:,end-t+1:end);
                    fprintf('using means for Bernoulli neurons...\n');
                    states = [sts1 sts2];
                case 'BP'
                    t = params.t;
                    sts1 = means(:,1:t);
                    fprintf('using means for Bernoulli neurons...\n');
                    sts2 = floor(means(:,(t+1):end));
                    fprintf('using modes for Poisson neurons...\n');
                    states = [sts1 sts2];
                case 'GB'
                    t = params.t;
                    sts1 = means(:,1:end-t);
                    fprintf('using means for Gaussian neurons...\n');
                    sts2 = round(means(:,end-t+1:end));
                    fprintf('using modes for Bernoulli neurons...\n');
                    states = [sts1 sts2];
                case 'GP'
                    t = params.t;
                    sts1 = means(:,1:t);
                    fprintf('using means for Gaussian neurons...\n');
                    sts2 = floor(means(:,(t+1):end));
                    fprintf('using modes for Poisson neurons...\n');
                    states = [sts1 sts2];
                case 'Bernoulli'
                    states = round(means);
                    if VERBOSE
                        fprintf('Using modes for Bernoulli neurons...\n');
                    end
                otherwise
                    states = means;
                    if VERBOSE
                        fprintf('using means for %s neurons...\n',FXN);
                    end
            end
    end
    
    if (layer == length(wts)/2) && DISP
        popDisplay(states,batch,params);
    end
end
Rout = states;


err = mean(mean((Rin - Rout).^2));

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function popDisplay(data,batch,params)

for i = 1:size(data,1)
    T = displayshape(data(i,:),params);
    PPCplot(cat(2,T{:}),params,'deepest-layer responses');
    drawnow;
    fprintf('batch %i, ex. %i\n', batch, i);
    pause()
end

end
%-------------------------------------------------------------------------%