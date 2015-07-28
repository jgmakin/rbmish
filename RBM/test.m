function ErrorStats = test(wts,params,varargin)
% test a trained RBM model

%-------------------------------------------------------------------------%
% Revised: 05/03/11
%   -pulled out estimatorStats-like stuff into a fxn estStatsCore.m
%   -pulled out NN-training/testing code core into fxn nnCore.m
% Revised: 05/03/11
%   -functionized
%   -created an option to train a new NN and test or just to test
%   -added some other cool stuff
% Created: 05/02/11
%   by JGM
%-------------------------------------------------------------------------%

% init
datagenargs = {};
ppcinfoargs = {};
updownargs = 'means';
iD = 1;
iP = 1;
% iU = 1;

if isempty(wts{1})
    error('the weights are empty; you probably forgot to assign them');
end

for i = 1:2:length(varargin)
    switch varargin{i}
        case {'visbias','propbias','prior'}
            datagenargs{iD} = varargin{i};
            datagenargs{iD+1} = varargin{i+1};
            iD=iD+2;
            
            ppcinfoargs{iP} = varargin{i};
            ppcinfoargs{iP+1} = varargin{i+1};
            iP=iP+2;
            
        case 'propagation'
            updownargs = varargin{i+1};
        case 'numsamples'
            params.smpls = varargin{i+1};
        otherwise
            datagenargs{iD} = varargin{i};
            datagenargs{iD+1} = varargin{i+1};
            iD=iD+2;
    end
end

% generate data
[pool,HADBEENCLOSED] = parallelInit;

[D0,S0] = DATAGENPP(1000,params,datagenargs{:});

% tensor -> matrix
[Di,Si] = longdata(D0,S0);
clear D0 S0;

% this is hacky...
% [Di xi] = killzeroinputs(Di,xi,params);

% up-down pass 
[avgErr,Do] = updown(Di,wts,params,updownargs);

% get statistics using normal (CoM) decoding
[statsL,statsN,ShatL] = estStatsCorePP(Si,params,'CoM',Di,Do);

% compute all the neutral-space error stats, and plot
%%%%
% xi is in the wrong format for PPCinfo....
%%%%
Shati = ShatL(:,:,:,1); % i.e., just the "in" estimates
ErrorStats = PPCinfo(Di,Shati,statsL,statsN,params,ppcinfoargs{:});
if HADBEENCLOSED, delete(pool); end

end













