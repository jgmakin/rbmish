% EFH   Exponential Family Harmonium training
%
%   EFH trains an exponential-family harmonium.  It requires the following
%   to have been defined outside this script:
%
%   EACHBATCHISATRAJ -- Boolean for dynamical data
%   RESTART     -- set to 1 if learning starts from beginning
%   Ntest       -- frequency in epochs of test-error updates
%   Npretrain   -- number of epochs of static training for dynamical data
%   Nbatches    -- number of batches/epoch, i.e. w/fixed learning params
%   Ncases      -- number of cases/batch, i.e. per weight change
%   NEFHs       -- how many EFHs in the stack of the deep belief network
%   iEFH        -- which in the stack of the deep belief network
%   getLatents  -- fxn handle for generating *data* latents
%   getData     -- fxn handle for generating observed data from latents
%
%   params      -- the master parameter structure (see params.m):
%       [[Might be good to try to get rid of all of these....]]
%       .numsUnits          -- units in each layer of the DBN
%       .typeUnits          -- *types* of units in each layer
%       .datatype           -- which data? (for initializing weights)
%       .mw, .mvb, .mhb
%       .bw, .bvb, .bhb
%       .kw, .kvb, .khb
%       .Ncdsteps
%       .algorithm          -- EFH, rEFH, TRBM, or RTRBM

%-------------------------------------------------------------------------%
% Revised: 12/08/16
%   -rewrote big chunks as part of Grand Revision.  The biggest change is
%   that, conceptually, the recurrent models no longer include recurrent
%   activity as part of their "data."  Practically, this means that EFH.m
%   concatenates on these activities, rather than overwriting the zeros (or
%   whatever) formerly placed there by generateData.m.
% Revised: 11/30/15
%   -cleaned up, rationalized, especially the plotting
% Revised: 01/15/15
%   -pulled out testEFHdecoding as a separate function
% Revised: 09/18/14
%   -renamed rbm.m -> EFH.m
% Revised: 05/14/13
%   -cleaned up
% Revised: 05/10/13
%   -added options for recurrency
% Revised: 09/14/10
%   -fixed tester to work with DBNs rather than just RBMs
% Revised: 05/25/10
%   -made number of contrastive-divergence steps variable
%   -incorporated rbmhidlinear
% Revised: 05/24/10
%   -fixed formatting
% modified by JGM
%-------------------------------------------------------------------------%



%---------------------------------- INIT ---------------------------------%
if RESTART == 1
    
    % extract/set params
    hidDstrbs = params.typeUnits{iEFH+1};
    visDstrbs = params.typeUnits{iEFH};
    hidNums = params.numsUnits{iEFH+1};
    visNums = params.numsUnits{iEFH};
    if strcmp(dataclass,'gpuArray')
        visNums = gpuArray(visNums);
        hidNums = gpuArray(hidNums);
    end
    [vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
        reinitializeEFH(iEFH,params.numsUnits,params.typeUnits,wts,...
        params.datatype,dataclass);
    RESTART=0; epoch=1;
    
    
    % rcrnt0 only *needs* to be set for recurrent models....
    if EACHBATCHISATRAJ
        rcrnt0 = zeros(size(hidbiases),'like',hidbiases);
    else
        rcrnt0 = zeros([Ncases,length(hidbiases)],'like',hidbiases);
    end
    
    algorithm = 'EFH';                  % these will be overwritten
    if (Npretrain>0)&&(iEFH==NEFHs)     %   after Npretrain.
        Ncdsteps = 5;                   %
    else
        Ncdsteps = params.Ncdsteps;
    end

    %%% inelegant: get rid of this by making params.mw (etc.) be a function
    %%% of iEFH as well as epoch
    if iEFH > 1
        params = setLearningSchedules(50,50,50,'hyperbolic',params,50,50);
        % NepochsMax = 40;
    end
    mwfunc = params.mw; mvbfunc = params.mvb;   mhbfunc = params.mhb;
    bwfunc = params.bw; bvbfunc = params.bvb;   bhbfunc = params.bhb;
    kwfunc = params.kw; kvbfunc = params.kvb;   khbfunc = params.khb;
    
    
    % for displaying
    Recons.DISP             = false(NepochsMax,1);
    Hiddens.DISP            = false(NepochsMax,1);
    Weights.DISP            = false(NepochsMax,1);
    ReconError.DISP         = true(NepochsMax,1);
    TestError.DISP          = false(NepochsMax,1);
    TestError.DISP(Ntest:Ntest:NepochsMax)=true;
    WeightNorm.DISP         = true(NepochsMax,1);
    WeightVelocityNorm.DISP = false(NepochsMax,1);
    
    
    if iEFH==1
        figmap = containers.Map({'Recons','Hiddens','Weights','ReconError',...
            'TestError','WeightNorm','WeightVelocityNorm'},{Recons,Hiddens,...
            Weights,ReconError,TestError,WeightNorm,WeightVelocityNorm});
    else % in case NepochsMax has changed...
         figmap('Recons') = setfield(figmap('Recons'),'DISP',Recons.DISP);
         figmap('Hiddens') = setfield(figmap('Hiddens'),'DISP',Hiddens.DISP);
         figmap('Weights') = setfield(figmap('Weights'),'DISP',Weights.DISP);
         figmap('ReconError') = setfield(figmap('ReconError'),'DISP',ReconError.DISP);
         figmap('TestError') = setfield(figmap('TestError'),'DISP',TestError.DISP);
         figmap('WeightNorm') = setfield(figmap('WeightNorm'),'DISP',WeightNorm.DISP);
         figmap('WeightVelocityNorm') = setfield(figmap('WeightVelocityNorm'),'DISP',WeightVelocityNorm.DISP);
    end
    figmap = EFHdisp(figmap,iEFH,0,wts,params);
   
end
%-------------------------------------------------------------------------%





%---------------------------------- LOOP ---------------------------------%
% cycle through training data
for epoch = epoch:NepochsMax
    errsum = 0;
    mw = mwfunc(epoch); mvb = mvbfunc(epoch); mhb = mhbfunc(epoch);
    bw = bwfunc(epoch); bvb = bvbfunc(epoch); bhb = bhbfunc(epoch);
    kw = kwfunc(epoch); kvb = kvbfunc(epoch); khb = khbfunc(epoch);
    
    % make or assemble new data every epoch
    batchdata = generateData(Nbatches*Ncases,getLatents,getData,dataclass);
    for iLayer = 1:(iEFH-1)
        batchdata = invParamMap(batchdata,wts{iLayer}(1:end-1,:),...
            wts{iLayer}(end,:),params.typeUnits{iLayer+1},...
            params.numsUnits{iLayer+1},params);
    end
    if EACHBATCHISATRAJ
        batchdata = shortdata(Nbatches,3,batchdata);
        batchdata = permute(batchdata,[3,2,1]);     % T x Nvis x Ntraj
    else
        batchdata = shortdata(Ncases,3,batchdata);  % Ntraj x Nvis x T
    end
    
    
    % loop through batches
    for iBatch = 1:Nbatches
        %%%fprintf('.');
        
        % if you're done w/pretraining, at the *top* EFH, and recurrent...
        %%% consider replacing Npretrain w/something based on norm(vishid)
        if (iEFH==NEFHs)&&(((epoch-1)*Nbatches + iBatch) == (Npretrain+1))
            
            % ...then the recurrent algorithms can kick in
            [visDstrbs,visNums,inputDstrbs,inputNums,inputInd,...
                vishid,vishidinc,visbiases,visbiasinc,algorithm,Ncdsteps,...
                mwfunc,mvbfunc,bwfunc,bvbfunc,kwfunc,kvbfunc,...
                mw,mvb,bw,bvb,kw,kvb] = recurrentize(visDstrbs,visNums,...
                hidDstrbs,hidNums,vishid,vishidinc,visbiases,visbiasinc,epoch,params);
        end
        
        
        % "positive" phase
        pvisstates = batchdata(:,:,iBatch);
        if EACHBATCHISATRAJ&&~strcmp(algorithm,'EFH')
            phidmeans = RNNforwardpass(pvisstates,rcrnt0,vishid,...
                hidbiases,hidDstrbs,hidNums,params);
            rcrntmeans = cat(1,rcrnt0,phidmeans(1:(end-1),:));
            pvisstates = cat(2,rcrntmeans,pvisstates);
        else
            if ~strcmp(algorithm,'EFH')
                if iBatch==1
                    pvisstates = cat(2,rcrnt0,pvisstates);
                else
                    pvisstates = cat(2,phidmeans,pvisstates);
                    % pvisstates = cat(2,phidstates,pvisstates);
                end
            end
            phidmeans = invParamMap(pvisstates,vishid,hidbiases,...
                hidDstrbs,hidNums,params);
        end
        phidstates = sampleT(phidmeans,hidDstrbs,hidNums,params);
        
        % negative phase
        switch algorithm
            case {'TRBM','RTRBM'}
                bh = hidbiases +...
                    pvisstates(:,1:(inputInd-1))*vishid(1:(inputInd-1),:);
                [qvisstates, qhidstates] = CDstepper(phidstates,...
                    vishid(inputInd:end,:),visbiases(inputInd:end),bh,...
                    hidDstrbs,inputDstrbs,hidNums,inputNums,Ncdsteps,params);
                qvisstates = cat(2,pvisstates(:,1:(inputInd-1)),qvisstates);
            otherwise
                [qvisstates, qhidstates] = CDstepper(phidstates,...
                    vishid,visbiases,hidbiases,...
                    hidDstrbs,visDstrbs,hidNums,visNums,Ncdsteps,params);
        end
        
        % the (negative of the) gradient "signal"
        p_dHdW = pvisstates'*phidstates;    q_dHdW = qvisstates'*qhidstates;    
        p_dHdhb = sum(phidstates);          q_dHdhb = sum(qhidstates);
        p_dHdvb = sum(pvisstates)';         q_dHdvb = sum(qvisstates)';
        gradientsigW  = p_dHdW  - q_dHdW;
        gradientsigbh = p_dHdhb - q_dHdhb;
        gradientsigbv = p_dHdvb - q_dHdvb;
        
        % does anything else accrue to the gradient?
        switch algorithm
            case 'RTRBM'
                [BPTTgradsigW,BPTTgradsigbh] = BPTT(pvisstates,phidmeans,...
                    (phidstates-qhidstates)*vishid(1:(inputInd-1),:)',...
                    vishid(1:(inputInd-1),:)',hidDstrbs);
                gradientsigW = gradientsigW + BPTTgradsigW;
                gradientsigbh = gradientsigbh + BPTTgradsigbh;
            case 'rEFH-BPTT'
                [BPTTgradsigW,BPTTgradsigbh] = BPTT(pvisstates,phidmeans,...
                    phidstates*vishid(1:Nhid,:)' + visbiases(1:Nhid)',...
                    vishid(1:Nhid,:)',hidDstrbs);
                gradientsigW = gradientsigW + BPTTgradsigW;
                gradientsigbh = gradientsigbh + BPTTgradsigbh;
            case 'ECcoherencesXXX'
                if iEFH==1
                    Nvis = visNums(1);
                    gradientsigW(1:Nvis/2,:) =...
                        repmat(mean(gradientsigW(1:Nvis/2,:)),[Nvis/2,1]);
                    gradientsigW((Nvis/2+1):end,:) =...
                        repmat(mean(gradientsigW((Nvis/2+1):end,:)),[Nvis/2,1]);
                    visbiases(1:Nvis/2) = mean(visbiases(1:Nvis/2));
                    visbiases((Nvis/2+1):end) = mean(visbiases((Nvis/2+1):end));
                end
        end
        
        % force sparse hidden units?
        if sparsitycost > 0
            if (epoch==1)&&(iBatch==1)
                vishidest   = p_dHdW;
                hidest      = p_dHdhb;
                visest      = p_dHdvb;
            else
                vishidest   = estimRate*p_dHdW  + (1-estimRate)*vishidest;
                hidest      = estimRate*p_dHdhb + (1-estimRate)*hidest;
                visest      = estimRate*p_dHdvb + (1-estimRate)*visest;
            end
            gradientsigW  = gradientsigW +...
                sparsitycost*(phidTarget*visest - vishidest);
            gradientsigbh = gradientsigbh +...
                sparsitycost*(Ncases*phidTarget - hidest);
        end
        
        % update weight/biases and their velocities simultaneously
        [vishid,vishidinc] = secondOrderWeightUpdate(vishid,...
            vishidinc,gradientsigW/Ncases,mw,bw,kw,visNums,params);
        [hidbiases,hidbiasinc] = secondOrderWeightUpdate(hidbiases,...
            hidbiasinc,gradientsigbh/Ncases,mhb,bhb,khb,hidNums,params);
        [visbiases,visbiasinc] = secondOrderWeightUpdate(visbiases,...
            visbiasinc,gradientsigbv/Ncases,mvb,bvb,kvb,visNums,params);
        
        
        % reconstruction error
        errsum = errsum + sum((pvisstates(:) - qvisstates(:)).^2);
    end
    fprintf('\n');
    
    % display
    figmap = EFHdisp(figmap,iEFH,epoch,wts,params,...
        vishid,hidbiases,visbiases,vishidinc,pvisstates,phidstates,...
        qvisstates,errsum/numel(batchdata));
    
    save rbmwts vishid hidbiases visbiases params epoch
    
    
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
% Version 1.000
%
% Code provided by Geoff Hinton and Ruslan Salakhutdinov
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with a
% note saying that the original programs are available from our web page.
% The programs and documents are distributed without any warranty, express
% or implied.  As the programs were written for research purposes only,
% they have not been tested to the degree that would be advisable in any
% important application.  All use of these programs is entirely at the
% user's own risk.

% This program trains Restricted Boltzmann Machine in which visible,
% binary, stochastic pixels are connected to hidden, binary, stochastic
% feature detectors using symmetrically weighted connections. Learning is
% done with 1-step Contrastive Divergence. The program assumes that the
% following variables are set externally:
%
%   maxepoch    -- maximum number of epochs
%   Nhid        -- number of hidden units
%   batchdata   -- the data that is divided into batches
%   (Ncases Ndims Nbatches)
%   RESTART     -- set to 1 if learning starts from beginning
%-------------------------------------------------------------------------%