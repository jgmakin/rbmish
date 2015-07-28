function [NNmean,NNcov,net] = nnCore(ShatN,StrueL,MOD,params,varargin)


% NB: THIS FUNCTION EXPECTS LONGDATA

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -XhatN -> ShatN, xi -> Si, all corresponding changes.
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% init params
Nexamples = size(ShatN,1);

% get (neutral-space) versions of the stimulus estimates, from pop. MOD
shatN = ShatN(:,:,strcmp(params.mods,MOD));

% get the "neutral space" version of the actual stimuli
sN = StrueL(:,:,strcmp(params.mods,params.NS));

% produce a neural net
if nargin > 4                                   % use a trained NN
    testvec = 1:Nexamples;
    net = varargin{1};
else                                            % train a new NN on errors
    trainvec = 1:Nexamples/2;
    testvec = Nexamples/2+1:Nexamples;
    e1 = shatN - sN;
    net = nnFit(shatN(trainvec,:)',e1(trainvec,:)');
end
    
% compute the NN output for all the data
eout = sim(net,shatN')';

% calculate the residual error for all data
e2 = (shatN - eout) - sN;

% get the stats on the testdata
NNmean = mean(e2(testvec,:));
NNcov = cov(e2(testvec,:));

end