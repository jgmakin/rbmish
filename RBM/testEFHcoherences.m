function testError = testEFHcoherences(R,X,Q,wts,params)


% Ns
Nvis = params.numsUnits{1}(end);
if strcmp(params.typeUnits{1}(end),'GammaFixedScale')
    Y = exp(R(:,1:Nvis));
else
    Y = R(:,(Nvis/2+1):Nvis);
end

% posterior probabilities and corresponding hidden-units params
[Pizq,Thyq] = getEMMposteriorProbs(R,wts,params);

% error in estimating hidden category
%%% get S from X?
%X = getCatProbs(S,params.typeUnits{2}{1});
%[fracCorrectTotal,Xntrp] = hiddenCatError(X,Pizq,1);

% "extra" plots
% updateSecretPlot(fracCorrectTotal*100,1002);
title('Percent Correct');
EMMparams = EFparams2ErlangParams(Pizq,Thyq);
plotEMMfit(Y,EMMparams,100,1051);
%updateSecretPlot(Xntrp,1003); title('Cross Entropy');
%MSE = Xntrp;                            % it has to be something
%%%%
foo = sampleT(Thyq,params.typeUnits{1}(end),params.numsUnits{1}(end),params);
MSE = mean((foo(:,end) - R(:,end)).^2);
%%%%


% put back into std units
testError = (det(MSE))^(1/size(MSE,1));


end