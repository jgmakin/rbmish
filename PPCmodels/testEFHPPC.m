function [testError,R1,Shat,Info] = testEFHPPC(R0,X,Q,wts,params,varargin)
% testEFHPPC    Standard decoding from an EFH trained on PPCs
%
% USAGES:
%   testError = testEFHPPC(R0,X,Q,wts,params);
%
%   [testError,R1,Shat,Info] = testEFHPPC(R0,X,Q,wts,params);
%
% Compute a scalar decoding error, testError, from the "responses" R0,
% latent variables X, data structure Q, DBN weights wts, and the standard
% parameters structure params.  More precisely, the testError is a kind of
% standard deviation of the (multi-dimensional) MSE of the neutral-space
% modality, calculated by taking a determinant and sending to the 1/Ndims
% power.
%
% This function also does some other possibly useful calculations along the
% way, including "updating" the inputs to R1 (via updown), and computing
% the sufficient statistics for the stimuli and gains from the "responses,"
% Shat and ttSpks.  These can be requested as outputs.
%
%   R:          (Nexamples x Nmods*N^Ndims)
%   X:          (Nexamples x Nstates x Nltnts)
%
%   Shat:       (Nexamples x Ndims x Nmods)
%   ttlSpks:    (Nexamples x Nmods)

%-------------------------------------------------------------------------%
% Revised: 01/05/17
%   -part of the Grand Revision
% Cribbed: 01/02/17
%   -from testEFHDecoding
%   by JGM
%-------------------------------------------------------------------------%


% default for updown
propagation = defaulter('updownargs','means',varargin{:});

% up and down, then decode (center of mass), then compute errors
if isfield(Q,'T')
    R1 = updownRDBN(R0,wts,params,Q.T,propagation);
else
    R1 = updownDBN(R0,wts,params,propagation);
end
[Shat,ttlSpks,err] = decodeDataPPC(R1,X,Q,params);

% compute MSE just on the errors of the "neutral-space" ind
errNS = err(:,:,strcmp(params.mods,params.NS));
MSE = cov(errNS) + mean(errNS)'*mean(errNS);

% put back into std units
testError = (det(MSE))^(1/size(MSE,1));

% if requested, get the posterior information matrix, too
if nargout > 3
   Info = GTPNposteriorInfo(ttlSpks,params); 
end

end