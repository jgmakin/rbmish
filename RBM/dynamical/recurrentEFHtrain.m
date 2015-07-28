% LearnDynamics
%   This script trains an "EFH filter" on data that encode a dynamical
%   system.  The dynamics themselves are set in setDynamics.m.  Then it
%   tests the trained network, with the classic test (error covariances
%   across all trials) as well as several visualization tools.

%-------------------------------------------------------------------------%
% Cribbed: 05/28/15
%   -from LearnDynamics.m
%   by JGM
%-------------------------------------------------------------------------%


% train an RBM-filter on dynamical data
clear; clc; 


load('dynamical/RTRBMtestdata150605.mat');
params = setParams;

%%% 1, 12, 20
Factor0 = 16;
Nfactors = 16;
Nruns = 1;
%%%

% testing data
[~,machine] = system('hostname');
params.machine = strtrim(machine);
%%% 
% params.dynamics.T = 1000;
% params.Ncases = 40;
%%%testData = getLDSdata(params);
params.Ncases = 400;
params.dynamics.T = 40000/params.Ncases;
%%% 
TESTDECODING = 1;
Ntest = 5;
EACHBATCHISATIMESTEP = 0;


% Ns/malloc
Nsensory = params.Nmods*params.N^params.Ndims;
errorStatMat = zeros(Nfactors,Nruns,'like',testData.Z);
ydataTensor = zeros(params.DBNmaxepoch/Ntest,Nfactors,Nruns);


for iFactor = Factor0:Nfactors
    

    params.numsUnits = [Nsensory + iFactor*Nsensory, iFactor*Nsensory];
    params.t = iFactor*Nsensory;
    params.numsUnits = [params.t+Nsensory, params.t];
    
    for iRun = 1:Nruns

        numsUnits = params.numsUnits;
        numRBMs = length(numsUnits)-1;
        paramDisplay(params);
        wts = cell(numRBMs*2,1);
        
        Nvis = params.numsUnits(1);
        datagenargs = {};
        
        % pretraining
        for i_rbm = 1:numRBMs
            Nhid = numsUnits(i_rbm+1);
            fprintf(1,'Pretraining Layer %i w/RBM: %d-%d \n',i_rbm,Nvis,Nhid);
            restart = 1;
            
            % train
            tic; rEFH; toc;
            
            % pack together weights for saving (hid => recog., vis => gener.)
            wts{i_rbm} = [vishid; hidbiases];
            wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases'];
            
            filename = 'EncoderWtsFile';
            save(filename,'i_rbm','numsUnits','wts','params','epoch');
            
            % for next time through
            Nvis = Nhid;
            %%% batchdata = batchphidmeans;
            
        end
        ydataTensor(:,iFactor,iRun) = gather(yvar); % ydata;
        
        % store the error variance
        [~,~,pEFH] = EFHfilter(testData,wts,params);
        pEFH.name = 'rEFH';
        
        % test
        [Ntraincases,~,Ntrainbatches] = size(batchdata);
        [Ntestcases,~,Ntestbatches] = size(testData.Z);
        params.dynamics.T = Ntestbatches;
        params.Ncases = Ntestcases;
        EFHstats = testDynamics(testData,params,0,pEFH);
        params.dynamics.T = Ntrainbatches;
        params.Ncases = Ntraincases;
        close all;
        errorStatMat(iFactor,iRun) = det(EFHstats(strcmp(params.mods,params.NS)).Cvrn);
    
        % store these weights if they're the best
        fprintf('iRun = %i, iFactor = %i\n',iRun,iFactor);
        currentStats = errorStatMat(Factor0:iFactor-1,:);
        currentStatVec = [currentStats(:); errorStatMat(iFactor,1:iRun)'];
        if errorStatMat(iFactor,iRun) == min(currentStatVec)
            bestwts = wts;
        end
       
    end

    % adjust the testing data for wts of the next size
%     Nhid = params.numsUnits(2);
%     [Ncases,Nvis,T] = size(D0);
%     inds = (Nhid+1):Nvis; %%% hard-coded for BP
%     PPCs = D0(:,inds,:);
%     fakedata = zeros(Ncases,Nhid+Nsensory,T);
%     D0 = cat(2,fakedata,PPCs);
   

end

save
