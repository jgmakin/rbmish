function varargout = marginalErrors(R,s,params,SILSCerrE,varargin)
% MARGINALERRORS    compute marginal error moments
%   MARGINALERRORS converts the conditional estimator moments into marginal
%   error moments.  Example:
%
%   [SINSMerrE SINSMerrT MINSMerrET MINSMerrT] =...
%       marginalErrors(xtrue,params,SILSCestE);
%
% Replaces local2common and computeOptimalCovs

%-------------------------------------------------------------------------%
% Revised: 08/12/11
%   -refined/fixed calculation of the trial-by-trial input variance.
% Revised: 06/01/11
%   -massive revision: functionized everything, cleaned up;
%   -fixed calculation on *input* covariances to use the law of total cov
% Revised: 04/12/11
%   -stuff
% Created: 04/08/11
%   by JGM
%-------------------------------------------------------------------------%


% init
Ndims = params.Ndims;
Nmods = params.Nmods; % 2
[Ncases,twom,Nbatches] = size(s);
n = 2; % length(varargin) + 1 %  nargin - 2;
bias = zeros(Ndims,Nmods);

% get biases, if there are any
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'visbias'
            LogInd = strcmp(params.mods,'Hand-Position');
            bias(:,LogInd) = varargin{i+1};
        case 'propbias'
            LogInd = strcmp(params.mods,'Joint-Angle');
            bias(:,LogInd) = varargin{i+1};
        case 'eyebias'
            LogInd = strcmp(params.mods,'Gaze-Angle');
            bias(:,LogInd) = varargin{i+1};
    end
end

% malloc
SINSCerrMu = zeros(n,Ndims,Nmods,Ncases*Nbatches);
MINSCerrMu = zeros(n,Ndims,Nmods,Ncases*Nbatches);
SINSCondCovSum = zeros(n,Ndims,Ndims,Nmods);
MINSCondCovSum = zeros(n,Ndims,Ndims,Nmods);

% computing tuning covariance
tuningCov = computetuningcovs(params);

% convert and average over the space
k = 0;
for iCase = 1:Ncases
    for iBatch = 1:Nbatches     
        
        k=k+1;
        sii = s(iCase,:,iBatch);
        
        SILSCerrT = PPCinputStats(R(iCase,:,iBatch),tuningCov,bias);
       
        [SINSCerrMu(:,:,:,k) SINSCondCov] = SICE(sii,params,...
            SILSCerrE,SILSCerrT);
        [MINSCerrMu(:,:,:,k) MINSCondCov] = MICE(sii,params,...
            SINSCerrMu(:,:,:,k),SINSCondCov);
        
        SINSCondCovSum = SINSCondCovSum + SINSCondCov;
        MINSCondCovSum = MINSCondCovSum + MINSCondCov;
    end
end

% convert conditional statistics to marginal statistics
SINSMerr = conds2marg(SINSCerrMu,SINSCondCovSum);
MINSMerr = conds2marg(MINSCerrMu,MINSCondCovSum);

% collect outputs
for j = 1:n 
    varargout{j} = SINSMerr(j,:);
    varargout{j+n} = MINSMerr(j,:);
end

% say
fprintf('\nUsing %d data...\n\n',Nbatches*Ncases);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [MINSCerrMu,MINSCondCov] = MICE(s,params,SINSCerrMu,SINSCondCov)
% multi-input conditional errors statistics

% init
Ndims = params.Ndims;
Nmods = params.Nmods; % 2;
n = size(SINSCondCov,1);
indices = reshape(1:Nmods*Ndims,Ndims,Nmods);
M = repmat(eye(Ndims),3,1);

% malloc
MINSCerrMu = zeros(n,Ndims,Nmods);
MINSCondCov = zeros(n,Ndims,Ndims,Nmods); 

% get prior
if isfield(params,'p0')
    prior.mu = params.p0.mu;
    prior.invcov = inv(params.p0.cov);
else
    prior.mu = zeros(Ndims,1);      % NB: dubious: assumes it'll never be used
    prior.invcov = zeros(Ndims);
end

% get multi-input error stats
for iMode = 1:Nmods
    for j = 1:n
        
        sL = s(indices(:,iMode));   
        %%%%%%%%%%%%%%%%%%%%%
        % shouldn't use the local one!!!
        %%%%%%%%%%%%%%%%%%%%%
        
        invCovs = [inv(squeeze(SINSCondCov(j,:,:,iMode)))...
            inv(squeeze(sum(SINSCondCov(j,:,:,:),4) - SINSCondCov(j,:,:,iMode)))...
            prior.invcov];
       
        % biases: pretty hacky and confusing---but right
        for iiMode = 1:Nmods
            b(:,iiMode) = (2*(iiMode<2)-1)*SINSCerrMu(j,:,iiMode)';
        end 
        B = (2*(iMode<2)-1)*[b(:,iMode); -(sum(b,2) - b(:,iMode))];
        biases = [B; prior.mu-sL'];
        %%%%%%%%%%%%%%%%%%%%%
        % but the prior is wrong in at least one modality
        %%%%%%%%%%%%%%%%%%%%%
        
        MINSCerrMu(j,:,iMode) = (invCovs*M)\invCovs*biases;        
        MINSCondCov(j,:,:,iMode) = inv(invCovs*M) -...      % SigmaPOST -...
            (invCovs*M)\prior.invcov/(invCovs*M);
        
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function MargErr = conds2marg(CondErrMu,CondCovSum)

% init
n = size(CondErrMu,1);
Ndims = size(CondErrMu,2);
Nmods = size(CondErrMu,3);
k = size(CondErrMu,4);


% malloc
MargErr = cell(n,Nmods);
for iMode = 1:Nmods
    for j = 1:n
        MargErr{j,iMode}.cov = zeros(Ndims);
        MargErr{j,iMode}.mu = zeros(Ndims,1);
    end
end

% loop
for iMode = 1:Nmods
    for j = 1:n
        CovOfMean = cov(squeeze(CondErrMu(j,:,iMode,:))');
        MeanOfCov = squeeze(CondCovSum(j,:,:,iMode))/k;
        MeanOfMean = mean(CondErrMu(j,:,iMode,:),4)';
        
        MargErr{j,iMode}.mu = MeanOfMean;
        MargErr{j,iMode}.cov = CovOfMean + MeanOfCov;
    end
end


end
%-------------------------------------------------------------------------%