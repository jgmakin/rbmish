function P = testPPCness(wts,params)
% Check to see if the hidden layer of one of your trained models is a PPC!
%
% It's not hard to check if a population is a PPC:  In the likelihood fxn
% of a PPC, the only interaction b/n the population (Y) and the stimulus 
% (X) comes through the inner product:
%
%       eta_x(x)'*T_x(y) = eta_x(x)'*y.  
%
% The same is therefore true of the posterior distribution, which will also
% [always? generally?] be an exponential family distribution, although its 
% sufficient statistics T_y(x) will not generally be linear.  Writing the 
% interaction term from the perspective of the posterior, we have:
%
%       eta_x(x)'*y = eta_y(y)'*T_y(x) 
%                  = (P'*y)'*T_y(x).
%       => eta_y(y) = P'y.
%
% Now, if you have arranged for the posterior to take a nice form---e.g., 
% Gaussian---then we can write out its natural parameters, eta_y(y)---e.g.,
% [Sigma^{-1}mu, -1/2*vect{Sigma^{-1}}].  Therefore, for each sample vector
% y, we can compute the corresponding posterior natural parameters, eta,
% and then regress to find the matrix P that relates them.
%
% If such a matrix does not exist, then the population Y is not a PPC.  For
% "input" populations that were designed to be PPCs, P obviously has to
% exist.  But for the hidden layer, it does not.  If we have already shown
% that the posterior conditioned on the hidden vector is (approximately)
% optimal, and this distribution is Gaussian, then PPCness of the hidden
% layer implies the existence of P---and contrapositively.
%
%
% USAGE:
%{
    load([getdir('data'),'RBMish/nonfinal/wts1DintegBBdropout92.mat']);
    testPPCness(wts,params)

    load([getdir('data'),'RBMish/EFHs/wts_2Dtwoarms_160226.mat']);
    testPPCness(wts,params)

    load([getdir('data'),'RBMish/EFHs/wts_1DrEFH_ManyXprmts.mat']);
    P = testPPCness(Allwts{1},params)
%}

%-------------------------------------------------------------------------%
% Created: 04/07/15
%   by JGM
%-------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%
% TO DO:
% (1) You need to make this (more) robust against redundancies in the
%   regression.  Maybe penalize it, or leave our more columns, or etc.
% (2) You need to account for the nonlinearity (kinematics)!!
%   Alternatively, you can test it on wts2Dtwoarms.mat.
%%%%%%%%%%%%%%%%%



%%%%%%%
% W = 100*W;
%%%%%%%

% train
[PPCs,trainPosteriors] = getPPCsAndPosteriors(wts,params);
[P,Rsq] = getLinearNPdecoders(PPCs,trainPosteriors);
cat(1,Rsq{:})

% test
[PPCs,testPosteriors,S] = getPPCsAndPosteriors(wts,params);
[~,Rsq,XpctHat] = getLinearNPdecoders(PPCs,testPosteriors,P{:});
fprintf('\n(Y1, Y2, decoder(Zbar), zbar) by natural parameters)\n');
cat(1,Rsq{:})


% check to see if this method performs as well as the standard one
checkErrorStats(S,XpctHat,testPosteriors.srcs,params);

    

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [PPCs,posteriors,S] = getPPCsAndPosteriors(wts,params)
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

if isfield(params,'dynamics')
    
    % make data
    T = 1000;
    Ntraj = 40;
    
    
    % get data
    [X,Q] = params.getLatents(Ntraj*T,dataclass);
    R0 = params.getData(X,Q);
    
    % decode inputs
    [Shat0,ttlSpks0] = decodeDataPPC(R0,X,Q,params);
    Info0 = GTPNposteriorInfo(ttlSpks0,params);
    pSENSORY = cumulantNeutralize(Shat0,Info0,params);
    
    % filter with rEFH
    [R1,Zparams] = updownRDBN(R0,wts,params,T);
    [Shat1,ttlSpks1] = decodeDataPPC(R1,X,Q,params);
    Info1 = GTPNposteriorInfo(ttlSpks1,params);
    pREFH = cumulantNeutralize(Shat1,Info1,params);
    
    % get the stimulus
    S = latents2stims(X,Q.latent2stim,params.mods,params.Ndims);
    S = reshape(S,[Ntraj,1,T]);

       
    %%%%%%%% under construction  
    load('dynamical\finalwts\LDSparamsEM2ndOrd1DrEFHManyXprmts.mat','Allparams');
    pEM2 = KFposteriorization(pSENSORY,Q,Allparams(1),params);
    pEM2.name = 'EM$^2$';
    %%%%
    % the posteriors---the targets of the regressions
    LDSparamsTrue = getLDSparams(params.dynamics);
    pOPT = KFposteriorization(pSENSORY,Q,LDSparamsTrue,params);
    pOPT.name = 'opt';
    
    
    %%%% how confident are the hidden units, really?
    % InfoTarg = (longdata(pOPT.Xpct)./longdata(pOPT.Cvrn))./longdata(pREFH.Xpct);
    %%%%
    posteriors.Xpct(:,:,1) = pOPT.Xpct;
    posteriors.Info(:,:,:,1) = 1./pSENSORY.Cvrn;
    %%% hard-coded for one-dimensional stimuli
    
    % LDSparamsREFH = fitLDStoKFobsvEsts(LDSdata,wts,params,'random');
    % LDSparamsREFH = fitLDStoKFobsvEsts(LDSdata,wts,params,LDSparamsREFH);
    load('dynamical\nonfinalwts\LDSparamsMirror.mat','LDSparamsREFH')
    pMirror = KFposteriorization(pSENSORY,Q,LDSparamsREFH,params);
    pMirror.name = 'mirror';
    posteriors.Xpct(:,:,2) = S; % pMirror.Xpct;
    pInfo2 = 1./pREFH.Cvrn;
    %%% hard-coded for one-dimensional stimuli
    posteriors.Info(:,:,:,2) = pInfo2;
    
    

    posteriors.Xpct(:,:,3) = pMirror.Xpct;
    posteriors.Info(:,:,:,3) = 1./pMirror.Cvrn; % pREFH.Cvrn;
    %%% hard-coded for one-dimensional stimuli
    posteriors.srcs = {pSENSORY.name,'hidsNonlin','hidsLinear'};
    
    
    % "PPCs," first and second layer---the inputs to the regressions
    PPCs{1} = R0;
    PPCs{2} = [pREFH.Xpct.*pInfo2,pInfo2];
    PPCs{3} = Zparams;
    %%%%%%%%
    
    
    %%%%
    borderRs = R0; 
    WRAPS = prod(borderRs(:,[1 15]),2); %%% not perfect...
    for i = 1:3
        PPCs{i} = PPCs{i}(~WRAPS,:);
    end
    % foo = 1./longdata(pOPT.Cvrn);
%     fakah(LDSparamsTrue.A,LDSparamsTrue.C,LDSparamsTrue.SigmaX,...
%         pSENSORY.Xpct,pSENSORY.Cvrn,pREFH.Xpct);
    
    posteriors.Xpct = posteriors.Xpct(~WRAPS,:,:);
    posteriors.Info = posteriors.Info(~WRAPS,:,:,:);
    S = S(~WRAPS,:);
    %%%%
    
else
    
    % generate training data (for the decoder), visibles and hiddens
    Nexamples = 20000;
    [Y,S] = generateData(Nexamples,params);
    Zparams = invParamMap(Y,wts{1}(1:end-1,:),wts{1}(end,:),...
        params.typeUnits{2},params.numsUnits{2},params);
    
    [unisensoryCumulants,multisensoryCumulants] = getPosteriorCumulants(Y,params);
    posteriors = assembleTestPosteriors(unisensoryCumulants,[],...
        multisensoryCumulants,[],[],params,'unisensory','optimal');
    
    % test the EFH against the optimal posterior
    posteriors.Xpct(:,:,end+1) = posteriors.Xpct(:,:,end);
    posteriors.Info(:,:,:,end+1) = posteriors.Info(:,:,:,end);
    posteriors.srcs{end+1} = 'EFH';
    %%%%
    posteriors.Xpct(:,:,1) = posteriors.Xpct(:,:,end);
    posteriors.Info(:,:,:,1) = posteriors.Info(:,:,:,end);
    posteriors.Xpct(:,:,2) = posteriors.Xpct(:,:,end);
    posteriors.Info(:,:,:,2) = posteriors.Info(:,:,:,end);
    %%%%
    
    
    
    % store the PPCs
    Nexamples = size(Y,1);
    Ysplit = reshape(Y,[Nexamples,params.N^params.Ndims,length(params.mods)]);
    
    for iPosterior = 1:length(posteriors.srcs)
        
        % the second one is in the wrong "space"!!  So you get Rsq < 1.0
        thisSrc = posteriors.srcs{iPosterior};
        switch thisSrc
            case {'Joint-Angle','Hand-Position'}
                PPCs{iPosterior} = Ysplit(:,:,strcmp(params.mods,thisSrc));
            case 'opt'
                PPCs{iPosterior} = Y;
            case 'EFH'
                PPCs{iPosterior} = Zparams;
            otherwise
                error('unsupported case --- jgm');
        end
    end
end
    

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [unisensoryCumulants,multisensoryCumulants] =...
    getPosteriorCumulants(Y,params)
    
% decode both populations as well as the optimal integrated estimate
[cntrOfMass, ttlSpks] = GTPNsuffstats(Y,params);
Info = GTPNposteriorInfo(ttlSpks,params);
%%%%%%
fprintf('\n\ndoing terrible clamping-at-edges thing to avoid complex numbers\n\n');
smin = shiftdim(params.smin,-1);
smax = shiftdim(params.smax,-1);
cntrOfMass = (cntrOfMass >= smin).*cntrOfMass + (cntrOfMass > smin).*smin;
cntrOfMass = (cntrOfMass <= smax).*cntrOfMass + (cntrOfMass > smax).*smax;
%%%%%%
unisensoryCumulants = cumulantNeutralize(cntrOfMass,Info,params);
multisensoryCumulants = gaussPosteriorization(unisensoryCumulants);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function testPosteriors = assembleTestPosteriors(unisensCmlnts0,...
    unisensCmlnts1,multisensCmlnts0,multisensCmlnts1,nnXpct,params,varargin)
% Assemble the posteriors to be tested and plotted!


% "audibles"
iPost = 0;
for iArg = 1:length(varargin)
    switch varargin{iArg}
        case 'unisensory'
            for iMod = 1:size(unisensCmlnts0.Xpct,3)
                iPost = iPost + 1;
                testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,iMod);
                testPosteriors.Info(:,:,:,iPost) = unisensCmlnts0.Info(:,:,:,iMod);
                testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{iMod};
            end
        case 'Joint-Angle'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,...
                strcmp(params.mods,'Joint-Angle'));
            testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{...
                strcmp(params.mods,'Joint-Angle')};
        case 'Hand-Position'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts0.Xpct(:,:,...
                strcmp(params.mods,'Hand-Position'));
            testPosteriors.srcs{iPost} = unisensCmlnts0.srcs{...
                strcmp(params.mods,'Hand-Position')};
        case 'optimal'
            iPost = iPost + 1;
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts0.Xpct;
            testPosteriors.Info(:,:,:,iPost) = multisensCmlnts0.Info;
            testPosteriors.srcs{iPost} = 'opt';
        case 'EFHuni'
            iPost = iPost + 1;
            indsN = strcmp(params.NS,unisensCmlnts1.srcs);
            testPosteriors.Xpct(:,:,iPost) = unisensCmlnts1.Xpct(:,:,indsN);
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHmulti'
            iPost = iPost + 1;
            fprintf('\ndecoding EFH by integrating "updated" populations\n');
            testPosteriors.Xpct(:,:,iPost) = multisensCmlnts1.Xpct;
            testPosteriors.srcs{iPost} = 'EFH';
        case 'EFHnn'
            iPost = iPost + 1;
            fprintf('\ndecoding EFH with ANN\n');
            testPosteriors.Xpct(:,:,iPost) = nnXpct;
            testPosteriors.srcs{iPost} = 'EFH';
        otherwise
            error('undefined posteriors! -- jgm\n');
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [P,Rsq,XpctHat] = getLinearNPdecoders(PPCs,posteriors,varargin)

for iPosterior = 1:length(posteriors.srcs)
    
    % the first one is in the wrong "space"!!  So you get Rsq < 1.0 
    thisPPC = PPCs{iPosterior};
    thisXpct = posteriors.Xpct(:,:,iPosterior);
    thisInfo = posteriors.Info(:,:,:,iPosterior);
    
    if isempty(varargin)
        [P{iPosterior},Rsq{iPosterior}] = fitLinearNPdecoder(thisPPC,...
            thisXpct,thisInfo);
        XpctHat = 0;
    else
        [P{iPosterior},Rsq{iPosterior},XpctHat{iPosterior}] =...
            fitLinearNPdecoder(thisPPC,thisXpct,thisInfo,varargin{iPosterior});
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [P,Rsq,XpctHat] = fitLinearNPdecoder(PPC,Xpct,Info,varargin)
% fit linear natural-parameter decoder


% the natural parameters of the posterior distribution, a Gaussian
eta1 = permute(tensorOp(permute(Info,[2,3,1]),permute(Xpct,[2,3,1])),[3,1,2]);
eta2 = -1/2*Info;

% PPC => eta = P'*y.  See labnotes ("GTPNs")
trialInds = sum(PPC,2)~=0;
unitInds = sum(PPC,1)~=0;
IN = [PPC(trialInds,unitInds) ones(sum(trialInds),1)];
OUT = [eta1(trialInds,:), eta2(trialInds,:)]; %%% works for matrices, too
%%%% IN = [PPC(trialInds,unitInds)./Info(trialInds,:), ones(sum(trialInds),1)];
%%%% OUT = Xpct;
%%%


if isempty(varargin)
    [Preduced,Rsq] = linregress(IN,OUT);
    P = zeros(size(PPC,2),size(OUT,2));
    
    for iOUT=1:size(OUT,2), P(unitInds,iOUT) = Preduced(1:end-1,iOUT); end
    P(end+1,:) = Preduced(end,:); % ye old column of ones
    
    XpctHat = 0;
else
    P = varargin{1};
    OUThat = IN*P([unitInds true],:);
    OUTbar = mean(OUT,1);
    Res = OUT - OUThat;
    SSerr = sum(Res.^2,1);
    SStot = sum((OUT - repmat(OUTbar,size(OUT,1),1)).^2);
    Rsq = 1 - SSerr./SStot;
    
    
    %%%%%%
    Nexamples = size(OUThat,1);
    CvrnHat = zeros(Nexamples,size(Info,2),size(Info,3));
    for iExample = 1:Nexamples
        InfoHat = reshape(-2*OUThat(iExample,(size(eta1,2)+1):end),...
            [size(Info,2),size(Info,3)]);
        CvrnHat(iExample,:,:) = inv(InfoHat);
    end
    XpctHat = permute(tensorOp(permute(CvrnHat,[2,3,1]),...
        permute(OUThat(:,1:size(eta1,2)),[2,3,1])),[3,1,2]);    
    %%% XpctHat = OUThat;
    %%%%%%
    
    
end




end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function checkErrorStats(S,XpctHat,srcs,params)

Stest = S(:,:,strcmp(params.mods,params.NS));
for iPosterior = 1:length(XpctHat)
    e = Stest - XpctHat{iPosterior};
    eStats.Xpct(:,iPosterior) = mean(e);
    eStats.Cvrn(:,:,iPosterior) = cov(e);
    eStats.tags(iPosterior).name = srcs{iPosterior};
    eStats.N(iPosterior) = sum(~isnan(e(:,1)));
end
dispErrStats(eStats,params.NS);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function S = unwrapStims(LDSdata,params)
% it's more like "use unwrapped versions" than "unwrap"

C = params.dynamics.C;
S = reshape(shortdata(size(LDSdata.S,1),3,longdata(LDSdata.Z)*C'),...
    size(LDSdata.S));

end
%-------------------------------------------------------------------------%



