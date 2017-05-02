function Rout = updownDBN(Rin,wts,params,varargin)
% UPDOWN   Propagation through a DBN
%   Rout = updownDBN(Rin,wts,params,varargin) propagates Rin through the DBN 
%   specified by wts and params, either (depending on your pov) all the way
%   through the "unfolded" network," or up and down the original one (tho' 
%   using different biases for each direction).
%
%   You can additionally specify the propagation method ('suffstats','MLE',
%   'modes','means'), where the default is something strange (probably a
%   combination of different propagation methods for different units).
%
%   Its original use in fine-tuning the deep autoencoder, but it turns out
%   to be useful in a lot of places.
%
%   The data must be of the for Nexamples x Nunits

%-------------------------------------------------------------------------%
% Revised: 01/04/17
%   -got rid of error as an output since it was only ever used to output
%   Rout.
%   -renamed updown.m -> updownDBN.m, in preparation for parallel with
%   updownRDBN.
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
Nlayers = length(wts)/2+1; % length(params.numsUnits);

propagation = 'means';
for i=1:length(varargin)
    switch varargin{i}
        case {'suffstats','MLE','modes','means','Nsamples'}
            propagation = varargin{i};
        case 'quiet'
            VERBOSE = 0;
    end
end
if VERBOSE
    fprintf('using method %s\n',propagation);
end


% cycle through layers (twice as many as the number of RBMs)
states = Rin;
for iLayer = 1:length(wts)
    
    % init
    jLayer = Nlayers-abs(Nlayers-iLayer-1); % really!
    dstrbs = params.typeUnits{jLayer};
    nums = params.numsUnits{jLayer};

    % the "means" ("moment parameterization")
    Mu = invParamMap(states,wts{iLayer}(1:end-1,:),wts{iLayer}(end,:),...
        dstrbs,nums,params);
%     if ~sum(strcmp(dstrbs,{'StandardNormal','Bernoulli','Poisson',...
%             'Binomial','Erlang','Gamma','Categorical','BernoulliDropout',}))
%         error('unrecognized NN fxn - jgm\n');
%     end
    
    switch propagation
        case 'suffstats'
            states = sampleT(Mu,dstrbs,nums,params);
        case 'means'
            for iGrp = 1:length(dstrbs)
                switch dstrbs{iGrp}
                    case 'BernoulliDropout'
                        states = invParamMap(states,20*wts{iLayer}(1:end-1,:),...
                            wts{iLayer}(end,:),dstrbs,nums,params);
                    otherwise
                        states = Mu;
                end
            end
        case 'modes'
            endinds = cumsum(nums);
            startinds = [1, endinds(1:end-1)+1];
            for iGrp = 1:length(dstrbs)
                inds = startinds(iGrp):endinds(iGrp);
                switch dstrbs{iGrp}
                    case 'StandardNormal'
                        states(:,inds) = Mu(:,inds);
                    case 'Poisson'
                        states(:,inds) = floor(Mu(:,inds));
                    case 'Bernoulli'
                        states(:,inds) = round(Mu(:,inds));
                    case 'Binomial'
                        n = params.Ntrials;
                        states(:,inds) = floor(Mu(:,inds)*(n+1)/n);
                    case 'Categorical'
                        keyboard
                        states(:,inds) = max(Mu(:,inds));
                    otherwise
                        error('modes don''t really make sense for this model')
                end
            end
        case 'Nsamples' % for "hidden layer decoding"
            %%% NB that this makes no sense for hidden units whose standard
            %%% parameters are not the mean!  E.g., Erlang, Gamma, etc.  In
            %%% fact, even for Poissons, this estimate is biased....
            if iLayer<Nlayers % only sample on the way up
                NNNN = params.smpls;
                states = sampleT(Mu,dstrbs,nums,params);
                for i=1:NNNN-1
                    states = states + sampleT(Mu,dstrbs,nums,params);
                end
                if VERBOSE
                    fprintf('using %i Bernoulli sample(s)...\n',NNNN);
                end
                states = states/NNNN;
            else
                states = Mu;
            end
        otherwise
            error('there is no longer a default for updown -- jgm');
%             switch dstrbs
%                 case 'Poisson'
%                     states = round(Mu); % means - 0.5;
%                     if VERBOSE
%                         fprintf('*rounding* Poisson means...\n');
%                     end
%                 case 'PB'
%                     t = params.t;
%                     sts1 = floor(Mu(:,1:end-t));
%                     fprintf('using modes for Poisson neurons...\n');
%                     sts2 = Mu(:,end-t+1:end);
%                     fprintf('using means for Bernoulli neurons...\n');
%                     states = [sts1 sts2];
%                 case 'BP'
%                     t = params.t;
%                     sts1 = Mu(:,1:t);
%                     fprintf('using means for Bernoulli neurons...\n');
%                     sts2 = floor(Mu(:,(t+1):end));
%                     fprintf('using modes for Poisson neurons...\n');
%                     states = [sts1 sts2];
%                 case 'GB'
%                     t = params.t;
%                     sts1 = Mu(:,1:end-t);
%                     fprintf('using means for Gaussian neurons...\n');
%                     sts2 = round(Mu(:,end-t+1:end));
%                     fprintf('using modes for Bernoulli neurons...\n');
%                     states = [sts1 sts2];
%                 case 'GP'
%                     t = params.t;
%                     sts1 = Mu(:,1:t);
%                     fprintf('using means for Gaussian neurons...\n');
%                     sts2 = floor(Mu(:,(t+1):end));
%                     fprintf('using modes for Poisson neurons...\n');
%                     states = [sts1 sts2];
%                 case 'Bernoulli'
%                     states = round(Mu);
%                     if VERBOSE
%                         fprintf('Using modes for Bernoulli neurons...\n');
%                     end
%                 otherwise
%                     fprintf('you didn''t really anticipate this case!!\n');
%                     keyboard
%                     states = Mu;
%                     if VERBOSE
%                         fprintf('using means for %s neurons...\n',dstrbs);
%                     end
%             end
    end
    
    if (iLayer == length(wts)/2) && DISP
        popDisplay(states,batch,params);
    end
end
Rout = states;


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