function [thhat,khat] = gammaDstrbML(logxbar,xbar)
% gammaDstrbML  Maximum-likelihood parameter estimation for the Gamma dstrb
%
% gammaDstrbML(logxbar,xbar) takes sufficient statistics mean(log(X)) and
% mean(X) from a data matrix X, each of whose rows are i.i.d samples from a
% Gamma distribution; and produces (a column vector of) estimates of the
% parameters, thhat and khat, for each of those distributions.

%-------------------------------------------------------------------------%
% Created: 03/24/16
%   by JGM
%-------------------------------------------------------------------------%


s = log(xbar) - logxbar;
khat = (3-s + sqrt((3-s).^2 + 24*s))./(12*s);

delkhat = -(log(khat) - psi(khat) - s)./(1./khat - psi(1,khat));
while max(abs(delkhat./khat)) > 0.0005
    khat = khat + delkhat;
    delkhat = -(log(khat) - psi(khat) - s)./(1./khat - psi(1,khat));
    fprintf('refining...\n');
end
fprintf('%f\n',max(abs(delkhat./khat)))
thhat = xbar./khat;


end
%-------------------------------------------------------------------------%