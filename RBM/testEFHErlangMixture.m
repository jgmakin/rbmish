function testError = testEFHErlangMixture(R,X,Q,wts,params)
% testErlangMixture
%
% USAGES:
%   [X,Q] = testEFHErlangMixture(Nexamples,yrclass,params)
%
%   [X,Q] = testEFHErlangMixture(Nexamples,yrclass,params,T)

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from testEFHDecoding.m
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nvis = params.numsUnits{1}(end);

% error in estimating Erlang parameters
X = getCatProbs(X,params.typeUnits{2}{1});
Thyp = X*[params.shapeparams,params.scaleparams];
[Pizq,Thyq] = getEMMposteriorProbs(R,wts,params);
err = Thyp - Thyq;
MSE = cov(err) + mean(err)'*mean(err);

% error in estimating hidden category
[fracCorrectTotal,Xntrp] = hiddenCatError(X,Pizq,1);

% "extra" plots
updateSecretPlot(fracCorrectTotal*100,1002); title('Percent Correct');
updateSecretPlot(Xntrp,1003); title('Cross Entropy');
%%%% You *could* also plot here the entropy of Z, a (low) floor
EMMparams = EFparams2ErlangParams(Pizq,Thyq);
Y = R(:,(Nvis/2+1):end);
plotEMMfit(Y,EMMparams,100,1051);
%updateSecretPlot(MSE(1,1),1003);
%updateSecretPlot(MSE(2,2),1003);

% put back into std units
testError = (det(MSE))^(1/size(MSE,1));

end