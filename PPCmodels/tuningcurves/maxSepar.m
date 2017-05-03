function [Vstar sv] = maxSepar(M)
% find the offset Vstar that maximizes the "separability" (dominance of the
% first singular value) of the matrix (M-Vstar), and return the first
% singular value of this matrix.
%
% This function was constructed to reproduce the analysis of Pena and
% Konisha 2001.

%-------------------------------------------------------------------------%
% Created: 08/24/12
%   by JGM
%-------------------------------------------------------------------------%

func = @(V)(sum(svd(M-V).^2) - max(svd(M-V).^2));
[Vstar,funcval,EXITFLAG] = fminbnd(func,-max(M(:)),max(M(:)));

if EXITFLAG==0
    fprintf('warning: fminbnd didn'' find a real sol''n -- jgm\n');
end

Mbar = M - Vstar;
s = svd(Mbar);
sv = s(1);

end