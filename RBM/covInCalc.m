function SINSMerr = covInCalc(R,S,params,varargin)
% find the error/estimator covariance of the different input modalities

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -now expects s0 (formerly x0) to be 4D rather than 3D.
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%


% init params
[Ncases,Ndims,Nmods,Nbatches] = size(S);

% malloc
SINSCerrMu = zeros(1,Ndims,Nmods,Ncases*Nbatches);
SINSCondCovSum = zeros(1,Ndims,Ndims,Nmods);

% compute tuning curve covariances
tuningCov = computetuningcovs(params);

% convert
tic
k = 0;
for iCase = 1:Ncases
    for iBatch = 1:Nbatches
        k=k+1;
        
        SILSCerrT = PPCinputStats(R(iCase,:,iBatch),tuningCov,zeros(Ndims,Nmods));   
        for iMod = 1:length(SILSCerrT)
            if isinf(SILSCerrT{iMod}.cov)
                fprintf('probably no spikes at (%i,%i)!\n',iCase,iBatch);
            end
        end
        [SINSCerrMu(:,:,:,k),SINSCondCov] =...
            SICE(squeeze(S(iCase,:,:,iBatch)),params,SILSCerrT);
        SINSCondCovSum = SINSCondCovSum + SINSCondCov;
    end
end
SINSMerr = conds2marg(SINSCerrMu,SINSCondCovSum);
toc

for iMod = 1:Nmods
    SINSMerr{iMod}.tags.src = 'single';
    SINSMerr{iMod}.tags.mod = params.mods{iMod};
    SINSMerr{iMod}.tags.space = 'neutral';
    SINSMerr{iMod}.tags.var = 'error';
    SINSMerr{iMod}.tags.epist = 'theoretical';
end

if nargin > 3
    if varargin{1} == 1
        dispErrCovs({SINSMerr},size(R,1)*size(R,3),params);
    end
end

end
%-------------------------------------------------------------------------%
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
for iMod = 1:Nmods
    for j = 1:n
        MargErr{j,iMod}.cov = zeros(Ndims);
        MargErr{j,iMod}.mu = zeros(Ndims,1);
    end
end

% loop
for iMod = 1:Nmods
    for j = 1:n
        CovOfMean = cov(squeeze(CondErrMu(j,:,iMod,:))');
        MeanOfCov = squeeze(CondCovSum(j,:,:,iMod))/k;
        MeanOfMean = mean(CondErrMu(j,:,iMod,:),4)';
        
        MargErr{j,iMod}.mu = MeanOfMean;
        MargErr{j,iMod}.cov = CovOfMean + MeanOfCov;
    end
end


end
%-------------------------------------------------------------------------%