tic
clear; clc;

if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

% params
Nxprmt  = 1;
Nruns   = 1;
Ntraj   = 40;
T       = 1000;

% construct the EM params
params          = setParams('datatype','LTI-PPC'); % ,'SensoryUnitType','Bernoulli');
Nstates         = size(params.dynamics.A,1);
NepochsMax      = params.NepochsMax;
datadstrb       = 'Normal';
muX0            = params.dynamics.muX0;
SigmaX0         = params.dynamics.SigmaX0;
muV0            = params.dynamics.muV0;
SigmaV0         = params.dynamics.SigmaV0;
mrgn            = 0.05; %%% bad hard-coding; see getLatentsLTI.m
%%% NB: for noninvertible C, this  can produce bad results!!
zmin            = pinv(params.dynamics.C)*params.smin(:);
zmax            = pinv(params.dynamics.C)*params.smax(:);
%%%
Ndims           = params.Ndims;
if strcmp(dataclass,'gpuArray')
    muX0 = gpuArray(muX0); muV0 = gpuArray(muV0);
end
getInitialState = @(Nparticles)(reshape(...
    [UniformNormalDiracSampler(muX0,SigmaX0,Ntraj*Nparticles,...
    zmin(1:Ndims),zmax(1:Ndims),mrgn),...
    UniformNormalDiracSampler(muV0,SigmaV0,Ntraj*Nparticles,[],[],mrgn)]',...
    [length(muX0)+length(muV0),Nparticles,Ntraj]));
getTrajs = @(yrclass)(EFHdata2LDSdata(40,@()(getLDStensor(yrclass))));





% the real model, or some reduced-order approximation?
Mstates = Nstates; % - 1;
getInitialStateWrapper = @(Nparticles)(reshape(eye(Mstates,Nstates)*...
    reshape(getInitialState(Nparticles),[Nstates,Nparticles*Ntraj]),...
    [Mstates,Nparticles,Ntraj]));



%%%
% bs = linspace(0,1,Nxprmt);
% ms = linspace(1,10,Nxprmt);
% ks = linspace(0,5,Nxprmt);
%%%

% malloc
allXNtrps = zeros(Nruns,Nxprmt);

for iXprmt = 1:Nxprmt
    XNtrpyBest = Inf;
    
    %%%
    % b=0.25; m=5; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -ks(iXprmt)/m*dt, -(b/m*dt-1)];
    % k=3; m=5; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -k/m*dt, -(bs(iXprmt)/m*dt-1)];
    % k=3; b=0.25; dt=0.05;
    % params.dynamics.A = [1.0000, dt; -k/ms(iXprmt)*dt, -(b/ms(iXprmt)*dt-1)];
    %%%
    
    for iRun = 1:Nruns
        
        HACKFLAG = 1;
        fprintf('experiment: %i, iteration: %i\n',iXprmt,iRun);        
        while HACKFLAG
            try
                LDSparamsEM = EM4LDS(Nstates-1,NepochsMax,datadstrb,...
                    getTrajs,'parameter initialization','random');
%                 LDSparamsEM = EM4LDS(Mstates,NepochsMax,datadstrb,...
%                     getTrajs,'parameter initialization','RandInit',...
%                     'inferencetype','particle',...
%                     'state initialization',getInitialStateWrapper);
                HACKFLAG = 0;
            catch ME
                keyboard
                fprintf('ugh, starting again\n');
            end
        end
        allXNtrps(iRun,iXprmt) = LDSparamsEM.XNtrp;
        allparams(iRun,iXprmt) = LDSparamsEM;
        
        if LDSparamsEM.XNtrp < XNtrpyBest
            XNtrpyBest = LDSparamsEM.XNtrp;
            LDSparamsEMthebest = LDSparamsEM;
        end
    end

    bestparams(iXprmt) = LDSparamsEMthebest;
end

toc


%%% save(['LDSparamsEM',params.datatype,params.dynamics.meta],'LDSparamsEM');




function [Y,X] = getLDStensor(yrclass,params)

[X,Q] = params.getLatents(40000,yrclass);
[R,Q] = params.getData(X,Q);
Y = decodeDataPPC(R,X,Q,params);

end