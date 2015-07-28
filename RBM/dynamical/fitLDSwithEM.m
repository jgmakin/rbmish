tic

% params
Nxprmt = 1;
Nruns = 20;
params = setParams;
params.dynamics.meta = 'RandInit';
%%% params.dynamics.meta = 'RandInitWithEC';
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;
params.Ncases = 40;
%%%
% bs = linspace(0,1,Nxprmt);
% ms = linspace(1,10,Nxprmt);
% ks = linspace(0,5,Nxprmt);
%%%

% malloc
allLLs = zeros(Nruns,Nxprmt);


for iXprmt = 1:Nxprmt
    LLbest = -Inf;
    
    %%%
    % b=0.25; m=5; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -ks(iXprmt)/m*dt, -(b/m*dt-1)];
    % k=3; m=5; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -k/m*dt, -(bs(iXprmt)/m*dt-1)];
    % k=3; b=0.25; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -k/ms(iXprmt)*dt, -(b/ms(iXprmt)*dt-1)];
    %%%
    
    for iRun = 1:Nruns
        
        fprintf('experiment: %i, iteration: %i\n',iXprmt,iRun);
        %%% LDSparamsEM = EM4LDS(size(params.dynamics.A,2),params);
        %%% LDSparamsEM = EM4LDS(size(params.dynamics.G,1),params);
        LDSparamsEM = EM4LDS(3,params);
        allLLs(iRun,iXprmt) = LDSparamsEM.LL;
        
        if LDSparamsEM.LL > LLbest
            LLbest = LDSparamsEM.LL;
            LDSparamsEMthebest = LDSparamsEM;
        end
    end

    Allparams(iXprmt) = LDSparamsEMthebest;
end

toc


%%% save(['LDSparamsEM',params.MODEL,params.dynamics.meta],'LDSparamsEM');