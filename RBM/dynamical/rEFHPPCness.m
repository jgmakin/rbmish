function rEFHPPCness

% how well can the posterior *natural parameters* be decoded linearly?
load('dynamical\finalwts\wts1DrEFHManyXprmts.mat')
wts = Allwts{1};
load('dynamical\nonfinalwts\LDSparamsMirror.mat','LDSparamsREFH')
load('dynamical\finalwts\LDSparamsEM2ndOrd1DrEFHManyXprmts.mat','Allparams')
LDSparamsEM = Allparams(1);
P = testPPCness(wts,params);



% testing and training data
[ZbarTest,Ztest,Stest,pEMtest,pMIRRORtest,pREFHtest] =...
    getXY(LDSparamsEM,LDSparamsREFH,wts,params,20);
[ZbarTrain,Ztrain,Strain,pEMtrain,pMIRRORtrain,pREFHtrain] =...
    getXY(LDSparamsEM,LDSparamsREFH,wts,params,40);



% only *plot* 1000 out of the 40000 data
inds = randperm(size(ZbarTest,1)*size(ZbarTest,3),1000);

% indices for random second-order interaction terms
Nhid = size(Ztest,2);
fooA = (1:Nhid)';
fooB = randperm(Nhid,Nhid)';
intinds = [fooA,fooB];


% regressions using Zbar
fprintf('\nLinear decoders of the stimulus and posterior mean:\n');
linearDecoder(ZbarTrain,ZbarTest,...
    cat(2,Strain,pEMtrain.Xpct,pMIRRORtrain.Xpct,pREFHtrain.Xpct,...
    pEMtrain.Xpct,pMIRRORtrain.Xpct,pREFHtrain.Xpct),...
    cat(2,Stest,Stest,Stest,pREFHtest.Xpct,...
    pEMtest.Xpct,pMIRRORtest.Xpct,pREFHtest.Xpct),...
    1,inds,...
    {'Zbar','S','xpct_{em}','xpct_{mirror}','xpct_{refh}',...
    'xpct_{em}','xpct_{mirror}','xpct_{refh}'},...
    {'Zbar','S','S','S','S','xpct_{em}','xpct_{mirror}','xpct_{refh}'});
    
fprintf('\nLinear decoders of the mean-to-variance ratio:\n');
linearDecoder(ZbarTrain,ZbarTest,...
    cat(2,pEMtrain.Xpct./pEMtrain.Cvrn,pEMtrain.Xpct./pMIRRORtrain.Cvrn,...
    pMIRRORtrain.Xpct./pEMtrain.Cvrn,pMIRRORtrain.Xpct./pMIRRORtrain.Cvrn,...
    pREFHtrain.Xpct./pEMtrain.Cvrn,pREFHtrain.Xpct./pMIRRORtrain.Cvrn),...
    cat(2,pEMtest.Xpct./pEMtest.Cvrn,pEMtest.Xpct./pMIRRORtest.Cvrn,...
    pMIRRORtest.Xpct./pEMtest.Cvrn,pMIRRORtest.Xpct./pMIRRORtest.Cvrn,...
    pREFHtest.Xpct./pEMtest.Cvrn,pREFHtest.Xpct./pMIRRORtest.Cvrn),...
    1,inds,...
    {'Zbar','xpct_{em}/vrnc_{em}','xpct_{em}/vrnc_{mirror}',...
    'xpct_{mirror}/vrnc_{em}','xpct_{mirror}/vrnc_{mirror}',...
    'xpct_{refh}/vrnc_{em}','xpct_{refh}/vrnc_{mirror}'},...
    {'Zbar','xpct_{em}/vrnc_{em}','xpct_{em}/vrnc_{mirror}',...
    'xpct_{mirror}/vrnc_{em}','xpct_{mirror}/vrnc_{mirror}',...
    'xpct_{refh}/vrnc_{em}','xpct_{refh}/vrnc_{mirror}'});

fprintf('\nLinear decoders of the inverse posterior variance:\n');
linearDecoder(ZbarTrain,ZbarTest,...
    cat(2,-1/2./pEMtrain.Cvrn,-1/2./pMIRRORtrain.Cvrn,-1/2./pREFHtrain.Cvrn),...
    cat(2,-1/2./pEMtest.Cvrn,-1/2./pMIRRORtest.Cvrn,-1/2./pREFHtest.Cvrn),...
    1,inds,...
    {'Zbar','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'},...
    {'Zbar','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'});

fprintf('\nNonlinear decoders of the inverse posterior variance:\n');
linearDecoder(... 
    cat(2,ZbarTrain,ZbarTrain(:,intinds(:,1),:).*ZbarTrain(:,intinds(:,2),:)),...
    cat(2,ZbarTest,ZbarTest(:,intinds(:,1),:).*ZbarTest(:,intinds(:,2),:)),...
    cat(2,-1/2./pEMtrain.Cvrn,-1/2./pMIRRORtrain.Cvrn,-1/2./pREFHtrain.Cvrn),...
    cat(2,-1/2./pEMtest.Cvrn,-1/2./pMIRRORtest.Cvrn,-1/2./pREFHtest.Cvrn),...
    1,inds,...
    {'Zbar,Zbar*Zbar','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'},...
    {'Zbar,Zbar*Zbar','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'});





% regressions using Z
fprintf('\nLinear decoders of the stimulus and posterior mean:\n');
linearDecoder(Ztrain,Ztest,...
    cat(2,Strain,pEMtrain.Xpct,pMIRRORtrain.Xpct,pREFHtrain.Xpct,...
    pEMtrain.Xpct,pMIRRORtrain.Xpct,pREFHtrain.Xpct),...
    cat(2,Stest,Stest,Stest,pREFHtest.Xpct,...
    pEMtest.Xpct,pMIRRORtest.Xpct,pREFHtest.Xpct),...
    1,inds,...
    {'Z','S','xpct_{em}','xpct_{mirror}','xpct_{refh}',...
    'xpct_{em}','xpct_{mirror}','xpct_{refh}'},...
    {'Z','S','S','S','S','xpct_{em}','xpct_{mirror}','xpct_{refh}'});
    
fprintf('\nLinear decoders of the mean-to-variance ratio:\n');
linearDecoder(Ztrain,Ztest,...
    cat(2,pEMtrain.Xpct./pEMtrain.Cvrn,pEMtrain.Xpct./pMIRRORtrain.Cvrn,...
    pMIRRORtrain.Xpct./pEMtrain.Cvrn,pMIRRORtrain.Xpct./pMIRRORtrain.Cvrn,...
    pREFHtrain.Xpct./pEMtrain.Cvrn,pREFHtrain.Xpct./pMIRRORtrain.Cvrn),...
    cat(2,pEMtest.Xpct./pEMtest.Cvrn,pEMtest.Xpct./pMIRRORtest.Cvrn,...
    pMIRRORtest.Xpct./pEMtest.Cvrn,pMIRRORtest.Xpct./pMIRRORtest.Cvrn,...
    pREFHtest.Xpct./pEMtest.Cvrn,pREFHtest.Xpct./pMIRRORtest.Cvrn),...
    1,inds,...
    {'Z','xpct_{em}/vrnc_{em}','xpct_{em}/vrnc_{mirror}',...
    'xpct_{mirror}/vrnc_{em}','xpct_{mirror}/vrnc_{mirror}',...
    'xpct_{refh}/vrnc_{em}','xpct_{refh}/vrnc_{mirror}'},...
    {'Z','xpct_{em}/vrnc_{em}','xpct_{em}/vrnc_{mirror}',...
    'xpct_{mirror}/vrnc_{em}','xpct_{mirror}/vrnc_{mirror}',...
    'xpct_{refh}/vrnc_{em}','xpct_{refh}/vrnc_{mirror}'});

fprintf('\nLinear decoders of the inverse posterior variance:\n');
linearDecoder(Ztrain,Ztest,...
    cat(2,-1/2./pEMtrain.Cvrn,-1/2./pMIRRORtrain.Cvrn,-1/2./pREFHtrain.Cvrn),...
    cat(2,-1/2./pEMtest.Cvrn,-1/2./pMIRRORtest.Cvrn,-1/2./pREFHtest.Cvrn),...
    1,inds,...
    {'Z','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'},...
    {'Z','vrnc_{em}','vrnc_{mirror}','vrnc_{refh}'});




% regressions with Z and (random) interaction terms
[fooA,fooB] = ndgrid(1:Nhid,1:Nhid);
fooA = triu(fooA,1);
fooB = triu(fooB,1);
intinds = [fooA(logical(fooA)),fooB(logical(fooB))];
% intinds = [1,2];

keyboard

[Rsq,MSE] = LMS4PPC(Ztrain,Strain,Ztest,Stest,intinds)
[Rsq,MSE] = LMS4PPC(Ztrain,-1/2./pEMtrain.Cvrn,Ztest,-1/2./pEMtest.Cvrn,intinds)

[Rsq,MSE] = LMS4PPC(ZbarTrain,Strain,ZbarTest,Stest,intinds)
[Rsq,MSE] = LMS4PPC(ZbarTrain,-1/2./pEMtrain.Cvrn,ZbarTest,-1/2./pEMtest.Cvrn,intinds)



%%Ztrain(:,intinds(:,1),:).*Ztrain(:,intinds(:,2),:))

fprintf('\nNonlinear decoders of the inverse posterior variance:\n');
linearDecoder(cat(2,Ztrain,Ztrain(:,intinds(:,1),:).*Ztrain(:,int2,:)),...
    cat(2,Ztest,Ztest(:,intinds(:,1),:).*Ztest(:,int2,:)),...,
    -1/2./pEMtrain.Cvrn,-1/2./pEMtest.Cvrn,1,inds,...
    {'Z','Z','vrnc_{em}','vrnc_{em}'}); % 20);
linearDecoder(cat(2,Ztrain,Ztrain(:,intinds(:,1),:).*Ztrain(:,int2,:)),...
    cat(2,Ztest,Ztest(:,intinds(:,1),:).*Ztest(:,int2,:)),...
    -1/2./pMIRRORtrain.Cvrn,-1/2./pMIRRORtest.Cvrn,1,inds,...
    {'Z','Z','vrnc_{mirror}','vrnc_{mirror}'}); % 20);
linearDecoder(cat(2,Ztrain,Ztrain(:,intinds(:,1),:).*Ztrain(:,int2,:)),...
    cat(2,Ztest,Ztest(:,intinds(:,1),:).*Ztest(:,int2,:)),...
    -1/2./pREFHtrain.Cvrn,-1/2./pREFHtest.Cvrn,1,inds,...
    {'Z','Z','vrnc_{refh}','vrnc_{refh}'}); % 20);



% %
% yrfactor = 200;
%
%
% for n = 1:yrfactor
%     inds = ((n-1)*(T/yrfactor)+1):(n*T/yrfactor);
%     foo = 1./squeeze(pMirrorTrain.Cvrn(:,:,:,inds));
%     caca = InfohatTrain(:,:,:,inds);
%     thing(n) = mean((foo(:)-caca(:)).^2);
% end
% figure(1); clf; plot(thing)



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function linearDecoder(Xtrain,Xtest,Ytrain,Ytest,t0,inds,trainvars,testvars)


X = longdata(Xtrain(:,:,t0:end));
%X = [X,ones(size(X,1),1)];
Y = longdata(Ytrain(:,:,t0:end));

[beta, RsqCV] = linrgsLOOCV(X,Y);
X = longdata(Xtest(:,:,t0:end));
%X = [X,ones(size(X,1),1)];
Y = longdata(Ytest(:,:,t0:end));


figure(3); clf;
Yhat = X*beta;
SSerr = sum((Y - Yhat).^2);
Ybar = mean(Y,1);
SStot = sum(bsxfun(@minus,Y,Ybar).^2,1);
Rsq = 1 - SSerr./SStot;
plot(Y(inds,:),Yhat(inds,:),'.')
axis equal
axis tight
% text(-.9,.9,['R$^2$ = ' sprintf('%2.4d\n',Rsq)])

for iOut = 1:size(Yhat,2)
    trainStr = sprintf('%-37s',['(',trainvars{1},'->',trainvars{iOut+1},'):']);
    testStr = sprintf('%-37s',['(',testvars{1},'->',testvars{iOut+1},'):']);
    fprintf(['TRAIN: ',trainStr,'TEST: ',testStr,'Rsq: %2.4f, '],Rsq(iOut));
    fprintf('MSE: %2.3d\n',SSerr(iOut)/size(Y,1));
end

pause(0.2);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Zbar,Z,S,pEM,pMIRROR,pREFH] =...
    getXY(LDSparamsEM,LDSparamsREFH,wts,params,Ntraj)

% Ns
params.Ncases = Ntraj;
T = params.dynamics.T;

% get inputs (Z) and outputs (S,eta,etc.) for regression
LDSdata = getLDSdata(params);
[Zbar,~,pREFH] = EFHfilter(LDSdata,wts,params);
Z = sampler(Zbar,params.typeUnits{2},params);
pMIRROR = KF4PPC(LDSdata,LDSparamsREFH,'mirror');
pEM = KF4PPC(LDSdata,LDSparamsEM,'EM2');
pREFH.Cvrn = reshape(pREFH.Cvrn,[Ntraj,1,T]);
pMIRROR.Cvrn = reshape(pMIRROR.Cvrn,[Ntraj,1,T]);
pEM.Cvrn = reshape(pEM.Cvrn,[Ntraj,1,T]);
S = reshape(LDSdata.S,[Ntraj,1,T]);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Rsq,MSE] = LMS4PPC(Ztrain,Ytrain,Ztest,Ytest,intinds)


% init
TOPLOT = 1;
[~,Nhid,T] = size(Ztrain);
th = randn(Nhid + size(intinds,1),1)/10;
thdot = zeros(size(th));
Nepochs = 500;
err = zeros(T,1);
Rsq = zeros(Nepochs,1);


% learned parameters
params.Ts = 1;
% Nvars = size(Ztrain,2) + length(intinds);
m = 38*2; % 968;
b = m/2; % 0.98*m; % m/6;
k = 0; % m/2000;

% prepare figure
if TOPLOT
    figure(9833);  % clf;
    hold on;
    plotHandle = plot(NaN,NaN);
    hold off;
end


for iEpoch = 1:Nepochs
    
    for t = 1:T
        
        % update
        y = Ytrain(:,:,t);
        %%% the ugly way is about 15% faster
        % z = cat(2,Ztrain(:,:,t),Ztrain(:,intinds(:,1),t).*Ztrain(:,intinds(:,2),t));
        % gradSig = LMS(th,z,y);
        gradSig = LMS2(th(1:Nhid),th((Nhid+1):end),Ztrain(:,:,t),...
            Ztrain(:,intinds(:,1),t).*Ztrain(:,intinds(:,2),t),y);
        %%%
        err(t) = norm(gradSig);
        [th,thdot] = secondOrderWeightUpdate(th,thdot,gradSig,m,b,k,'none',params);
        if err(t)>100000, fprintf('\nerror is huge: %d\n',err(t)); break; end
        
        % print errors
        if 0
            set(0,'CurrentFigure',figure(9833));
            hold on;
            set(plotHandle,'XData',1:t,'YData',err(1:t),'color','r');
            hold off;
            title(num2str(iEpoch))
            drawnow
        end
    end


    [Ntraj,Nhid,T] = size(Ztest);
    SSerr = 0;
    for t = 1:T
        y = Ytest(:,:,t);
        % z = cat(2,Ztest(:,:,t),Ztest(:,intinds(:,1),t).*Ztest(:,intinds(:,2),t));
        % SSerr = SSerr + sum((y - z*th).^2);
        SSerr = SSerr + sum((y - Ztest(:,:,t)*th(1:Nhid) -...
            Ztest(:,intinds(:,1),t).*Ztest(:,intinds(:,2),t)*th((Nhid+1):end)).^2);
    end
    Ybar = mean(longdata(Ytest),1);
    SStot = sum((longdata(Ytest) - repmat(Ybar,Ntraj*T,1)).^2);
    Rsq(iEpoch) = 1 - SSerr/SStot;
    MSE = SSerr/Ntraj/T;

    if TOPLOT
        set(0,'CurrentFigure',figure(9833));
        hold on;
        set(plotHandle,'XData',1:iEpoch,'YData',Rsq(1:iEpoch),'color','r');
        hold off;
        % title(num2str(iEpoch))
        drawnow
    else
        fprintf('Epoch: %i; Rsq: %2.4f; MSE: %2.4d\n',iEpoch,Rsq(iEpoch),MSE);
    end
end

Rsq = Rsq(end);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function delTh = LMS(th,z,y)

delTh = z'*(y-z*th)/size(y,1);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function delTh = LMS2(th1,th2,Z1,Z2,y)

Z1y = Z1'*y;
Z2y = Z2'*y;
Zth = Z1*th1 + Z2*th2;
fooA = Z1y - Z1'*Zth;
fooB = Z2y - Z2'*Zth;
delTh = [fooA;fooB]/size(y,1);

end
%-------------------------------------------------------------------------%











