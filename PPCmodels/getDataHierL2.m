function [R,Q] = getDataHierL2(X0,Q0,wts0,params0,params)

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   by JGM
%   from generateData.m
%-------------------------------------------------------------------------%


%%% TO DO:
% (1) Work out how to save the gains created here (Q.G)

% init
visDstrbs = params.typeUnits{1};
visNums = params.numsUnits{1};
endinds = cumsum(visNums);
startinds = [1, endinds(1:end-1)+1];
Q.Q0 = Q0;
Q.params0 = params0;


% the lower-level EFH
hidDstrbs0  = params0.typeUnits{2};
hidNums0    = params0.numsUnits{2};
smpls       = params0.smpls;

% malloc
R = zeros(size(X0,1),sum(visNums),'like',X0);

% make data differently depending on whether or not its from below
for iGrp = 1:length(visDstrbs)
    switch visDstrbs{iGrp}
        case 'Poisson'

            % useful things
            Nexamples = size(X0,1);
            uniformSmplFxn = @(M,xmin,xmax)(UniformNormalDiracSampler(...
                zeros(length(xmin),1,'like',X0),Inf,M,xmin,xmax,0));

            % what is normally done in getLatents[]
            Q.G = uniformSmplFxn(Nexamples,params.gmin,params.gmax);
            Nmods = length(params.mods);
            Q.latent2stim = cell(1,Nmods);
            for iMod = 1:Nmods
                f = setSensoryTransformations(params.mods{iMod},...
                    params.NS,params0.roboparams);
                Q.latent2stim{iMod} = @(X)(f(X(:,:,1),X(:,:,2)));
            end
            params.typeUnits{1} = visDstrbs(iGrp);
            params.numsUnits{1} = visNums(iGrp);
            
            % now get PPC responses at this (upper) level
            R(:,startinds(iGrp):endinds(iGrp)) = getDataPPC(X0,Q,params);
            
            
        case 'Bernoulli'
            if isa(X0,'gpuArray')
                for iEFH = 1:length(wts0)
                    wts0{iEFH} = gpuArray(wts0{iEFH}); 
                end
            end
            
            Q.R1and2 = getDataPPC(X0,Q0,params0);
            means = invParamMap(Q.R1and2,wts0{1}(1:end-1,:),...
                wts0{1}(end,:),hidDstrbs0,hidNums0,params0);
            %%%
            states = sampleT(means,hidDstrbs0,hidNums0,params0);
            for i = 1:smpls-1
                states = states + sampleT(means,hidDstrbs0,hidNums0,params0);
            end
            %%% Or:
            % params0.Ntrials = smpls;
            % states = sampleT(means,{'Binomial'},hidNums0,params0);
            %%% but somehow, the for-loop is faster, at least with the
            %%% GPU--maybe because of memory issues
            R(:,startinds(iGrp):endinds(iGrp)) = states/smpls;
            
        otherwise, error('unexpected unit types in hierarchy -- jgm');
    end
end

end