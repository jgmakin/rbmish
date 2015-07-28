function [ErrorStats,net] = nnDecode(wts,MOD,params,varargin)
% NNDECODETRAIN     Neural network for decoding
%   NNDECODETRAIN learns a neural network which maps the center-of-mass
%   from the MOD population to its corresponding decoding error.  The
%   network is then validated on test data: 
%
%       Run the test data through the NN.  Then subtract the outputs, which 
%       correspond to errors, from those same test data; these are the
%       "corrected" responses.  Now subtract the actual RBM inputs from
%       these responses---this is the overall error.  Find its mean (the
%       bias) and covariance and plot on top of the optimal.
%

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S and all associated changes
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
netvararg = {};
datagenargs = {};
ppcinfoargs = {};
iD = 1;
iP = 1;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'pretrained'
            netvararg{1} = varargin{i+1};
        case {'visbias','propbias','eyebias'}
            datagenargs{iD} = varargin{i};
            datagenargs{iD+1} = varargin{i+1};
            iD=iD+2;
            
            ppcinfoargs{iP} = varargin{i};
            ppcinfoargs{iP+1} = varargin{i+1};
            iP=iP+2;
        otherwise
            datagenargs{iD} = varargin{i};
            datagenargs{iD+1} = varargin{i+1};
            iD=iD+2;
    end
end

% generate data
[D0,S0] = DATAGENPP(1000,params,datagenargs{:});

% up-down pass
[Di,Si] = longdata(D0,S0);
clear D0 x0;
[avgErr,Do] = updown(Di,wts,params,'means');

% get statistics using normal (CoM) decoding
[statsL,statsN, ~, ShatNo] = estStatsCorePP(Si,params,'CoM',Di,Do);
clear blank;

% use a neural net to get better estimators
[NNstats{3}.mu NNstats{3}.cov net] = nnCore(ShatNo,Si,MOD,params,netvararg{:});
NNstats{1} = NNstats{3};
%%%%%%%%%%%%%%%%%%
NNstats{2} = NNstats{3};
%%%%%%%%%%%%%%%%%%
% rather than re-do it for the other MOD


% compute all the neutral-space error stats, and plot
ppcinfoargs{iP} = 'nn';
ppcinfoargs{iP+1} = NNstats;
ErrorStats = PPCinfo(Di,xx,statsL,statsN,params,ppcinfoargs{:});

end