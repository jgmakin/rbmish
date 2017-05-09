function XNtrp = GaussianCrossEntropy(innov,S,mattype,VERBOSE)
% GaussianCrossEntropy  from the "innovation" and marginal covariance
%
% USAGE:
%   XNtrpZ = GaussianCrossEntropy(innovZ,SigmaZ)
%
% ...
%
% Used in FactorAnalysis.m and KalmanFilter.m.
%
% NB: innovZ should have size (Nz x Nsamples)

%-------------------------------------------------------------------------%
% Updated: 05/02/16
%   -changed to skip XNtrp calculation when S is not positive definite
% Cribbed: 04/22/16
%   from FactorAnalysis.m
%   by JGM
%-------------------------------------------------------------------------%

N = size(innov,1);
[Sonehalf,p] = chol(S);
if p > 0
    if VERBOSE
        fprintf('cholesky decompsition failed...\n');
        %%%fprintf('trying eigenvalue decomposition...\n');
    end
    %%%D = eig(S);
    %%%cvrncost = sum(log((real(D+50*eps))));
    switch mattype, case 'info',XNtrp = inf; case 'cvrn',XNtrp = -inf; end
else
    cvrncost = 2*sum(log(diag(Sonehalf)));
    switch mattype
        case 'info'
            XNtrp = (mean(sum((S*innov).*innov)) + N*log(2*pi) - cvrncost)/2;
        case 'cvrn'
            XNtrp = (mean(sum((S\innov).*innov)) + N*log(2*pi) + cvrncost)/2;
        otherwise
            error('the mattype must either be ''info'' or ''cvrn'' -- jgm');
    end    
end
if VERBOSE, fprintf('cross entropy: %.3f\n',XNtrp); end

end