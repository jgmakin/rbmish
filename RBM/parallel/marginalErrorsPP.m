function [SINSMerr,MINSMerr] = marginalErrorsPP(R,S,params,SILSCerrE,varargin)
% MARGINALERRORS    compute marginal error moments
%   MARGINALERRORS converts the conditional estimator moments into marginal
%   error moments.  Example:
%
%   [SINSMerrE SINSMerrT MINSMerrET MINSMerrT] =...
%       marginalErrors(xtrue,params,SILSCestE);
%
% Replaces local2common and computeOptimalCovs
%
% NB: This function expects to get longdata!!

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S and all associated changes
% Revised: 12/10/13
%   -fixed ancient "bug" (?): no need for longdata in here: this function 
%   expects to get longdata!  This now causes an error b/c you changed the
%   way longdata.m works.
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
[Nexamples,Ndims,Nmods] = size(S);
Nargs = 2; % length(varargin) + 1 %  nargin - 2;
bias = zeros(Ndims,Nmods);
p0 = [];

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
        case 'prior'
            p0 = varargin{i+1};
    end
end

% malloc
SINSCondCovSum = zeros(Ndims,Ndims,Nmods,Nargs);
MINSCondCovSum = zeros(Ndims,Ndims,Nmods,Nargs);

% computing tuning covariance
tuningCov = computetuningcovs(params);

% convert and average over the space
[pool,HADBEENCLOSED] = parallelInit;
for iExample = 1:Nexamples
    
    s = shiftdim(S(iExample,:,:),1);
    SILSCerrT = PPCinputStats(R(iExample,:),tuningCov,bias);
    
    [SINSCerrMu(iExample,:,:,:),SINSCondCov]=SICE(s,params,...
        SILSCerrE,SILSCerrT);
    [MINSCerrMu(iExample,:,:,:),MINSCondCov]=MICE(s,params,...
        SINSCerrMu(iExample,:,:,:),SINSCondCov,p0);
    
    SINSCondCovSum = SINSCondCovSum + SINSCondCov;
    MINSCondCovSum = MINSCondCovSum + MINSCondCov;
end
if HADBEENCLOSED, delete(pool); end

% convert conditional statistics to marginal statistics
SINSMerr = conds2marg(SINSCerrMu,SINSCondCovSum);
MINSMerr = conds2marg(MINSCerrMu,MINSCondCovSum);

% say
fprintf('\nUsing %d data...\n\n',Nexamples);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [MINSCerrMu,MINSCondCov] = MICE(shat,params,SINSCerrMu,SINSCondCov,p0)
% multi-input conditional errors statistics

% init
[nequals1,Ndims,Nmods,Nargs] = size(SINSCerrMu);
M = repmat(eye(Ndims),3,1);  
%%% for summing three square matrices stacked side by side

% malloc
MINSCerrMu = zeros(Ndims,Nmods,Nargs);
MINSCondCov = zeros(Ndims,Ndims,Nmods,Nargs); 


% malloc
E = NaN(Ndims,Nmods);

% get multi-input error stats
for iMod = 1:Nmods
    for iArg = 1:Nargs
        
        % get prior
        if isempty(p0)
            e0 = zeros(Ndims,1);
            W0 = zeros(Ndims);
        else 
            e0 = p0.mu - shat(:,iMod) + squeeze(SINSCerrMu(1,:,iMod,iArg))';
            % This is the actual "error" given by the prior "estimator": it
            % thinks s is at p0.mu, and s is actually at shat - e = s.  NB
            % that shat is *local* (i.e., s is in its own space for each
            % modality), b/c that's what SICE expects.  Therefore, since
            % the error is in the neural space, this computation breaks for
            % the non-neutral-space modality/ies.
            %
            % But then, the prior is only for one modality, anyway, and
            % it's hard to translate it into the other.  This whole
            % approach is really bad, then, and should (one day) be
            % changed....
            W0 = inv(p0.cov);  
        end
        
        
        W1 = inv(SINSCondCov(:,:,iMod,iArg));
        W2 = inv(squeeze(sum(SINSCondCov(:,:,:,iArg),3)) - SINSCondCov(:,:,iMod,iArg));    
        invCovs = [W1 W2 W0];
        
        % biases: pretty hacky and confusing---but right
        for jMod = 1:Nmods
            E(:,jMod) = (2*(jMod<2)-1)*SINSCerrMu(1,:,jMod,iArg)';
        end
        e = (2*(iMod<2)-1)*[E(:,iMod); -(sum(E,2) - E(:,iMod))];
        biases = [e; e0];
        
        % ``integrate'' the errors
        invCovSum = invCovs*M;
        MINSCerrMu(:,iMod,iArg) = invCovSum\invCovs*biases;        
        MINSCondCov(:,:,iMod,iArg) = inv(invCovSum) -...    % SigmaPOST -...
            invCovSum\W0/invCovSum;
        
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function MargErr = conds2marg(CondErrMu,CondCovSum)

% Ns
[Nexamples,Ndims,Nmods,Nargs] = size(CondErrMu);

% malloc
MargErr = cell(Nargs,Nmods);
for iMod = 1:Nmods
    for iArg = 1:Nargs
        MargErr{iMod,iArg}.cov = zeros(Ndims);
        MargErr{iMod,iArg}.mu = zeros(Ndims,1);
    end
end

% loop
for iMod = 1:Nmods
    for iArg = 1:Nargs
        CovOfMean = cov(CondErrMu(:,:,iMod,iArg));
        MeanOfCov = CondCovSum(:,:,iMod,iArg)/Nexamples;
        MeanOfMean = mean(CondErrMu(:,:,iMod,iArg),1);
        
        MargErr{iMod,iArg}.mu = MeanOfMean(:);
        MargErr{iMod,iArg}.cov = CovOfMean + MeanOfCov;
    end
end


end
%-------------------------------------------------------------------------%