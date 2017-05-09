% EFHrestart
% restarts rbm training from last saved data (rbmwts)

%-------------------------------------------------------------------------%
% Revised: 09/18/14
%   -renamed rbmrestart.m -> EFHrestart.m
% Revised: 09/15/14
%   -again cleaned up, tested
% Revised: 05/14/13
%   -"cleaned up," but didn't test
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% clear; clc;
% load rbmwts;
% load('results\nonfinalwts\MCDwts140915.mat')

%%%%% THINGS YOU SHOULD EDIT
params.NepochsMax = 400;
NepochsMax = params.NepochsMax;
params = setLearningSchedules(100,100,100,'hyperbolic',params,1); 
%%%%%%


% THINGS DBN.m WOULD HAVE DONE
iRBM = 1; %%% by assumption
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
numsUnits = params.numsUnits;
numRBMs = length(numsUnits)-1;
Ncases = params.Ncases;
Nbatches = params.Nbatches;
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
Nvis = sum(visNums);
Nhid = sum(hidNums);
datagenargs = {}; %%% depends on what model you're loading
TESTDECODING = 1;
DISPLAYTESTS = 0;
Ntest = 5;


% THINGS EFH.m WOULD HAVE DONE (if restart == 1)
% extract/set params
datagenargs = {'dbndepth',iRBM,'dbnwts',wts};
[vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(iRBM,params.numsUnits,params.typeUnits,wts,...
    params.datatype,dataclass);
Ncdsteps = params.Ncdsteps;
Recons.DISP = false(NepochsMax,1);
Hiddens.DISP = false(NepochsMax,1);
Weights.DISP = false(NepochsMax,1);
ReconError.DISP = true(NepochsMax,1);
TestError.DISP=false(NepochsMax,1);TestError.DISP(Ntest:Ntest:NepochsMax)=true;
WeightNorm.DISP = true(NepochsMax,1);
WeightVelocityNorm.DISP = false(NepochsMax,1);

if iRBM==1
    figmap = containers.Map({'Recons','Hiddens','Weights','ReconError',...
        'TestError','WeightNorm','WeightVelocityNorm'},{Recons,Hiddens,...
        Weights,ReconError,TestError,WeightNorm,WeightVelocityNorm});
else
    %%%%% until you fix these---altho' at the moment these don't do
    %%%%% anything....
    Weights.DISP=false(NepochsMax,1);
    Recons.DISP=false(NepochsMax,1);
end
figmap = EFHdisp(figmap,datagenargs,iRBM,0,wts,params);




% NOW OVERWRITE WITH THE LOADED WEIGHTS
vishid = wts{iRBM}(1:end-1,:);
hidbiases = wts{iRBM}(end,:);
visbiases = wts{numRBMs*2-iRBM+1}(end,:)';

% run 
EFH

% store the resulting weights!
wts{iRBM} = [vishid; hidbiases];
wts{numRBMs*2-iRBM+1} = [vishid'; visbiases'];