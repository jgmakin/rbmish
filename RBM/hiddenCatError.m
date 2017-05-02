function [fracCorrectTotal,Xntrp] = hiddenCatError(Cats,Pizq,TOPLOT)
% Hidden-category errors using max assignment
%
% USAGE:
%   [Pizq,Thyq] = getEMMposteriorProbs(wts,params,testData.R);
%   fracCorrectTotal = hiddenCatError(longdata(testData.S),Pizq);

%-------------------------------------------------------------------------%
% Cribbed: 06/02/16
%   from testEFHdecoding.m
%   by JGM
%-------------------------------------------------------------------------%


% assign to categories based on biggest probability
Cathats = catAssign(Cats,Pizq);

% error wrt true categories [[this is not a great measure...]]
fracCorrectTotal = mean(sum(round(Cathats) == Cats,2)==size(Cathats,2));
ErrorRateByCat = mean(Cathats - Cats);
fprintf('Percent error by category: ');
for iCat = 1:length(ErrorRateByCat)
    fprintf('%.2f, ',ErrorRateByCat(iCat)*100);
end
fprintf('\n\n');

% cross entropy (a better measure)
Xntrp = -mean(log(Cathats(logical(Cats))));

% plot?
if TOPLOT
    figure(1001)
    imagesc([Pizq(1:200,:),Cats(1:200,:)])
end

end
