function MargErrorStats = PPCinfo(R,Shat,statsL,statsN,params,varargin)
% PPCINFO Compute and display useful information-theoretic quantities
%   PPCINFO takes the generating stimulus th and the stats structure
%   produced by estimatorStats.m---along with the master parameters---and
%   produces for display:
%
%       -the MLE error covariances for vis and prop
%       -the actal error covariances of the input and output populations
%       -the determinants of these variarances
%       -(mutual information)
%       -plots of the covariances
%
% NB: The input s is used *only* to compute the (manipulator) Jacobian.
% That means that it shouldn't really matter whether you use the
% uncorrelated or correlated vis and th (and indeed in practice it doesn't
% appear to matter).

%-------------------------------------------------------------------------%
% Revised: 06/17/14
%   -added case for varargin 'prior'
% Revised: 08/12/11
%   -corrected calculation of theoretical optimum in marginalErrors, which
%   no longer computes the theoretical input variances via Fisher info.
%   Hence PPCexpectedFI was cut, which in turn allowed the elimination of
%   all (exlicit) references to the gain
%   -this also require making the data matrix an argument
% Revised: 05/12/11
%   -got rid of varargout; now returns *all* the marginal error statistics
% Revised: 05/05/11
%   -changed plotting functions to accomodate an extra input, the neural
%   net's statistics
% Revised: 04/08/11
%   -lots of changes
%       --replaced local2common and computeOptimalCovs by marginalErrors.m.
%        This corrects some mistakes in calculating the optimal covariance
%        from the input covs (see labnotes.pdf)
%       --added MAP estimate stuff
%       --added an input for the "neutral-space" stats, b/c it's not
%       obvious how to convert the RBM's output marginal error covariance
%       from prop space to vis space (whereas it's straightforward to
%       calculated it if the RBM's *estimates* were converted before
%       computing stats)
% Revised: 02/11/11
%   -added plotting for 1D variances
%   -added varargin for bias
% Revised: 01/28/11
%   -added varargouts
% Revised: 01/27/11
%   -radically reworked in the service of computing 3-input covariances
%   correctly
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme
% Revised: 11/29/10
%   -changed to accomodate x in *true* coordinates
% Revised: 11/10/10
%   -changed PPCCovMLE to operate on a tuning covariance matrix rather than
%   a mere variance
% Revised: ??/??/10
%   -changed a lot of stuff!  Imp.: implemented true, general Jacobian
%   calculation---only to replace it again w/the specific, fast version
% Created: 08/23/10
%   -by JGM
%-------------------------------------------------------------------------%


%--------------------------------- INIT ----------------------------------%
% init params
Nmods = params.Nmods;
Ndims = params.Ndims;
SILSerrE = statsL(:,1);     % single-input, local-space, marg/cond errors
MINSMerrE = statsN(:,2);    % multi-input, neutral-space, marginal errors
bias = zeros(Ndims,Nmods);
MargErrorStats = cell(1,4);

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'visbias'
            LogInd = strcmp(params.mods,'Hand-Position');
            bias(:,LogInd) = varargin{i+1};
            SILSerrE{LogInd}.mu = SILSerrE{LogInd}.mu + bias(:,LogInd);
        case 'propbias'
            LogInd = strcmp(params.mods,'Joint-Angle');
            bias(:,LogInd) = varargin{i+1};
            SILSerrE{LogInd}.mu = SILSerrE{LogInd}.mu + bias(:,LogInd);
        case 'eyebias'
            LogInd = strcmp(params.mods,'Gaze-Angle');
            bias(:,LogInd) = varargin{i+1};
            SILSerrE{LogInd}.mu = SILSerrE{LogInd}.mu + bias(:,LogInd);
        case 'prior'
            fprintf('testing with a nonflat prior...\n');
            %%% anything else??
        case 'nn'
            nnStats = varargin{i+1};
            MargErrorStats{5} = nnStats;         %%% hacky!
    end
end
%-------------------------------------------------------------------------%


%------------------- COMPUTE MARGINAL ERROR COVARIANCES ------------------%
[SINSMerr,MINSMerr] = marginalErrorsPP(R,Shat,params,SILSerrE,varargin{:});

for iMod = 1:Nmods                                         % *** see (3) ***
    MINSMerrE{iMod}.mu = MINSMerrE{iMod}.mu + SINSMerr{iMod,1}.mu;
    if length(MargErrorStats) == 5
        MargErrorStats{5}{iMod}.mu =...
            MargErrorStats{5}{iMod}.mu'+SINSMerr{iMod,1}.mu;
    end
end
%-------------------------------------------------------------------------%


%-------------------------------- DISPLAY --------------------------------%
% display error covariances
fprintf('single-input error covariances\n');
printErrCovs(SINSMerr(:,2),'theoretical');
printErrCovs(SINSMerr(:,1),'empirical');

fprintf('multi-input error covariances\n');
printErrCovs(MINSMerr(:,2),'theoretical I (pure)');
printErrCovs(MINSMerr(:,1),'theoretical II (mixed)');
printErrCovs(MINSMerrE,'empirical');

if Ndims > 1
    % display error determinants
    fprintf('single-input error determinants\n');
    fprintf('   theoretical:        empirical:\n')
    for iMod = 1:size(SINSMerr,1)
        fprintf('   %d         %d\n',...
            det(SINSMerr{iMod,2}.cov),det(SINSMerr{iMod,1}.cov));
    end
    fprintf('multi-input error determinants\n');
    fprintf('   theoretical I:    theoretical II:      empirical\n')
    for iMod = 1:size(MINSMerr,1)
        fprintf('   %d      %d      %d\n',det(MINSMerr{iMod,2}.cov),...
            det(MINSMerr{iMod,1}.cov),det(MINSMerrE{iMod}.cov));
    end
end

% plot ellipses
MargErrorStats{1} = SINSMerr(:,1);
MargErrorStats{2} = SINSMerr(:,2);
MargErrorStats{3} = MINSMerrE;
MargErrorStats{4} = MINSMerr(:,1);
% MargErrorStats{5} is NNstats, if available
MargErrorStats{end+1} = MINSMerr(:,2);


% add some identifying tags
for iMod = 1:Nmods
    for iStat = 1:length(MargErrorStats)
        switch iStat
            case {1,2}
                MargErrorStats{iStat}{iMod}.tags.src = 'single';
                MargErrorStats{iStat}{iMod}.tags.name = params.mods{iMod};
                MargErrorStats{iStat}{iMod}.tags.mod = params.mods{iMod};
            case {3,4}
                MargErrorStats{iStat}{iMod}.tags.src = 'multiple';
                MargErrorStats{iStat}{iMod}.tags.name = 'opt';
                MargErrorStats{iStat}{iMod}.tags.mod = 'opt';
            otherwise
                MargErrorStats{iStat}{iMod}.tags.src = 'multiple';
                MargErrorStats{iStat}{iMod}.tags.name = 'EFH';
                MargErrorStats{iStat}{iMod}.tags.mod = 'multisensory';
        end
        
        MargErrorStats{iStat}{iMod}.tags.space = 'neutral';
        MargErrorStats{iStat}{iMod}.tags.var = 'error';
        if (iStat == 1) || (iStat == 3)
            MargErrorStats{iStat}{iMod}.tags.epist = 'empirical';
        else
            MargErrorStats{iStat}{iMod}.tags.epist = 'theoretical';
        end
    end
end

% display the stats of the inputs (1), outputs (3), and theor. opt. (3)
dispErrCovs(MargErrorStats([1,3,end]),size(R,1)*size(R,3),params);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function printErrCovs(cellstruct,str)

fprintf(['  ',str,'\n'])
for i = 1:length(cellstruct)
    disp(cellstruct{i}.cov);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function H = gaussEntr(S,p)
% the entropy of a normal distribution with covariance S

H = 1/2*log((2*pi*exp(1))^p*det(S));

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function I = PPCMI(cov,range,params)
% optimal decoding

% init
Ndims = params.Ndims;                                   % dimensionality
r = range(:,2) - range(:,1);
respLength0 = [params.respLength params.respLength];
respLength = scalefxn(respLength0,[0 0],respLength0,[0 0],r);
respA = prod(respLength);

% stimulus entropy
Hth = -log(1/respA);

% entropy of normal distr.
Hgauss = gaussEntr(cov,Ndims);

% mutual information (H(th) - E(Hvgauss))
% *** see (2) *** %
I = Hth - Hgauss;                               % >= actual I(x,xhat)

end
%-------------------------------------------------------------------------%


% *** (2) ***
% the "equivocation" is the average (over stimuli) information lost by
% Gaussian fluctations around the mean.  The real losses will be no larger
% than this (b/c the normal distr. has max. entropy for a given mean and
% variance).
% Becase H[Gauss] is not a fxn of theta, it is its own expected value.

% *** see (3) ***
% Admittedly confusing: You compute the empirical errors in estError w.r.t.
% possibly biased "true" stimuli, so here you have to add on the input
% biases.  That is, say e.g. your mechanism for implementing a discrepancy
% b/n modalities was to add b--say a positive scalar/vector---to prop.  The
% error-computing fxn estError doesn't know this, and will find (if it's
% working correctly) a *negative* bias.  To plot this with vis at [0;0],
% you have to add on the prop input bias.  Yeah.



%-------------------------------------------------------------------------%
% pos = scalefxn(s(i,1:2,j),patchmin,patchmax,vmin,vmax);
% rsq = pos(1)^2 + pos(2)^2;
% z1 = (L1^2 + L2^2 - rsq)/(2*L1*L2);
% z2 = (rsq + L1^2 - L2^2)/(2*L1*sqrt(rsq));
%
% mInvJ = [(1-z1^2)^(-1/2)*pos(1)/(L1*L2),...
%     (1-z1^2)^(-1/2)*pos(2)/(L1*L2);...
%     pos(2)/rsq + (1-z2^2)^(-1/2)*(4*pos(1)*sqrt(rsq) - (rsq + L1^2 + L2^2)*rsq^(-1/2)*L1*2*pos(1)),...
%     -pos(1)/rsq + (1-z2^2)^(-1/2)*(4*pos(2)*sqrt(rsq) - (rsq + L1^2 + L2^2)*rsq^(-1/2)*L1*2*pos(2))];
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% compute mutual-information upper bounds
% Ivi = PPCMI(stats.covVi,[vmin,vmax],params);      % >= actual I(v,vhat_i)
% Ipi = PPCMI(stats.covPi,[pmax,pmin],params);      % >= actual I(p,phat_i)
% Ivo = PPCMI(stats.covVo,[vmin,vmax],params);      % >= actual I(v,vhat_o)
% Ipo = PPCMI(stats.covPo,[pmax,pmin],params);      % >= actual I(p,phat_o)
% IF = PPCMI(covMLEvv,[vmin,vmax],params);      % >= Ivi,Ipi,etc; (opt.)
% Imax = PPCMI(covOpt,params);
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% display mutual-information bounds
% fprintf('information upper bounds:\n');
% fprintf('I(v,vhat_i): %f\nI(p,phat_i): %f\nI(v,vhat_o): %f\n',Ivi,Ipi,Ivo);
% fprintf('I(p,phat_o): %f\nI_fisher:    %f\nImax:
% %f\n',Ipo,IF,Imax);
%-------------------------------------------------------------------------%