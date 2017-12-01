function [testError,XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,...
    Bv,LDSparamsObs] = testEFHBMI(Rtest,Xtest,Qtest,wts,params,varargin)
% testEFHBMI    Decoding from an EFH trained on BMI data
%
% USAGE:
%   testError = testEFHBMI(Rtest,Xtest,Qtest,wts,params);
%
%   [testError,XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,...
%       Bv,LDSparamsObs] = testEFHBMI(Rtest,Xtest,Qtest,wts,params);
%
%   [testError,XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,...
%       Bv,LDSparamsObs] = testEFHBMI(Rtest,Xtest,Qtest,wts,params,...
%       Rtrain,Xtrain,Qtrain,SStot);

%-------------------------------------------------------------------------%
% Cribbed:
%   01/02/17
%   from testEFHDecoding
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) Consider computing MSE in two dimensions, a la the PPC models (so one
% for each of pos, vel, and acc), rather than R^2....

TOPLOT = 0;

% the decoder needs some training data, too
if isempty(varargin)
    [Xtrain,Qtrain] = params.getLatents([],class(wts{1}),...
        'sequencelength','singlesequence');
    [Rtrain,Qtrain] = params.getData(Xtrain,Qtrain);
    SStot = sum((Xtest - mean(Xtest)).^2);
else
    [Rtrain,Xtrain,Qtrain,SStot] = deal(varargin{:});
end

% observed variables
Vtrain = LDSobsfxn(Rtrain,Qtrain,wts,params);
clear Rtrain Qtrain
Vtest = LDSobsfxn(Rtest,Qtest,wts,params);
clear Rtest Qtest

% decode
[XhatStatic,RsqStatic,XhatDynamic,RsqDynamic,Bv,LDSparamsObs] =...
    filterDecoder(Vtrain,Xtrain,Vtest,Xtest,1,SStot,params.Nmsperbin);

% report -<SNR>
testError = mean(10*log10(1 - RsqStatic(:))); 

% plot?
if TOPLOT
    figure(1); clf;
    subplot(1,2,1); imagesc(wts{1}(1:params.numsUnits{1}(1),:)');
    colorbar
    subplot(1,2,2); imagesc(wts{1}((params.numsUnits{1}(1)+1):end-1,:)');
    colorbar
end
        
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function obsvs = LDSobsfxn(R0,Q,wts,params)

% Ns
Nunits = params.numsUnits{end};
Nlayers = length(params.numsUnits);
Mlayers = length(wts)/2+1;
Nsamples = size(R0,1);
T = Q.T;
clear Q;

% up and down through the rEFH
[R1,Z0] = updownRDBN(R0,wts,params,T);
if Nlayers == Mlayers
    % for the final layer of the rEFH, use all units...
    U1 = invParamMap(Z0,wts{Nlayers}(1:end-1,1:sum(Nunits)),...
        wts{Nlayers}(end,1:sum(Nunits)),params.typeUnits{end},...
        Nunits,params);
    obsvs = [R1,U1,Z0,ones(Nsamples,1,'like',R0)];
else
    % ...for the early layers, don't use the hidden units
    obsvs = R1;
end


end
%-------------------------------------------------------------------------%