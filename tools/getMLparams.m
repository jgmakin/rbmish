function MLparams = getMLparams(t,MLparams)
% The maximum likelihood parameters, written in terms of the (minimal?)
% sufficient statistics, t.  
%
% NB that these statistics are the *expected outer products," not the
% covariance matrices.

%-------------------------------------------------------------------------%
% Revised: 08/28/14
%   -eliminated code for "noisily observed controls" (this is taken care of
%   elsewhere)
%   -removed flags for 'muX' and  'muY', replaced with checks for the size
%   of the sufficient statistics.
%   -added calculation for matrix D in Y = CX + DU + noise
% Revised: 01/15/14
%   -accommodated (noisily observed) control inputs
% Revised: 11/26/13
%   -changed to assume that the sufficient stats included a column of ones
%   -changed again to make the presence of muX and muY signal "biased"
%   (nonzero) transition and emission noise (respectively).  This required
%   adding MLparams (LDSparams) as an input argument.
% Revised: 11/25/13
%   -got rid of M = (M+M')/2, which seemed gratuitous.
% Cribbed: 10/30/13
%   -from getLDSparams
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%

% FLAGS

% Ns
[Nstates,Nxp] = size(t.XfXp); 
Nxy = size(t.YX,2);

% regression coefficients
betaXY = t.YX/t.XX;
betaXX = t.XfXp/t.XpXp;

% A, B, and muX
MLparams.A = betaXX(:,1:Nstates);
if Nxp > Nstates
    MLparams.muX = betaXX(:,end);
    if Nxp > (Nstates+1)
        MLparams.B = betaXX(:,(Nstates+1):(end-1));
    end
end
MLparams.SigmaX = t.XfXf - betaXX*t.XfXp';

% C, D, and muY
MLparams.C = betaXY(:,1:Nstates);
if Nxy > Nstates
    MLparams.muY = betaXY(:,end);
    if Nxy > (Nstates+1)
        MLparams.D = betaXY(:,(Nstates+1):(end-1));
    end
end
MLparams.SigmaY = t.YY - betaXY*t.YX';


% prior parameters
MLparams.mu0 = t.mu0;
MLparams.Info0 = inv(t.x0x0 - t.mu0*t.mu0');
% --- (1) --- %

end


% --- (1) ---% 
% A little confusing, but this is right.
%
% For symmetry with the other terms, learnfullyobservedLDS.m and EM4LDS.m
% return an expected outer product, E[X0*X0'] (called t.x0x0), rather than
% the covariance of the initial state, Cov[X0], which is the required term
% (or its inverse, Info0).  These functions also return E[X0], so that one
% can calculate the initial-state covariance via:
%
%   Cov_x[X0] = E_x[X0*X0'] - E_x[X0]*E_x[X0'].
%
% In learnfullyobservedLDS.m, the terms on the left are calculated directly
% (since X0 is observed), and is the covariance *across trajectories*.  In
% EM4LDS.m, they are computed via the law of total expectation:
%
%   E_x[X0] = E_y[E_{x|y}[X0|Y]];
%
%   E_x[X0*X0] = E_y[E_{x|y}[X0*X0'|Y]]
%              = E_y[ Cov_{x|y}[X0|Y] + E_{x|y}[X0|Y]*E_{x|y}[X0'|Y] ]
%              = E_y[Cov_{x|y}[X0|Y]] + E_y[E_{x|y}[X0|Y]*E_{x|y}[X0'|Y]].
%
% The same applies, mutatis mutandis, to Cov[U].
%
% Ghahramani and Hinton (1996) provide a different equation (Eq'n 25):
%
%   Cov_x[X0] = P0 - E_x[X0]*E_x[X0'] + Cov_y[E_{x|y}[X0|Y]].
%
% They had defined P0 := E_{x|y}[X0*X0'|y]--but this presumably they meant
% this for the one-trajectory case, where "averaging" out Y from this eq'n
% amounts to doing nothing.  But if they meant this P0 to signify the
% average, i.e. P0 ?= E_y[E_{x|y}[X0*X0'|y]] = E_x[X0*X0'], then this eq'n
% is wrong, since it is just our first eq'n, *plus an extra covariance
% term*.  How now?





