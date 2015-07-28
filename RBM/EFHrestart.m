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


% THINGS setParams.m WOULD HAVE DONE
path(path,'../ee125/fxns');
path(path,'parallel');
path(path,'results');
path(path,'dynamical');
path(path,'tuningcurves');
path(path,'retired');
path(path,'scratch');
path('../utils',path);  % contains your tex.m
path(path,'../tools')



% THINGS DBN.m WOULD HAVE DONE
i_rbm = 1; %%% by assumption
[~,machine] = system('hostname');
params.machine = strtrim(machine);
numsUnits = params.numsUnits;
numRBMs = length(numsUnits)-1;
Ncases = params.Ncases;
Nbatches = 1000;
Nvis = numsUnits(1);
Nhid = numsUnits(i_rbm+1);
datagenargs = {}; %%% depends on what model you're loading
TESTDECODING = 1;
DISPLAYTESTS = 0;
Ntest = 5;


% THINGS EFH.m WOULD HAVE DONE (if restart == 1)
% extract/set params
HIDFXN = params.typeUnits{i_rbm+1};
VISFXN = params.typeUnits{i_rbm};
% maxepoch = params.DBNmaxepoch;
% amass = params.amass;


%%%%% THINGS YOU SHOULD EDIT
epoch = 51;
% amass = 1.0;
maxepoch = 90;
%%%%%



datagenargs = [datagenargs,{'dbndepth',i_rbm,'dbnwts',wts}];
[vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(i_rbm,params.numsUnits,wts);
batchposhidmeans = zeros(Ncases,Nhid,Nbatches,'like',params.mw);
allErrors = zeros(maxepoch,1,'like',params.mw);
restart=0; erravg=0; tErravg=inf; trErravg=inf; counter=0;

% plot errors
if TESTDECODING
    
    if isfield(params,'dynamics')
        testData = getLDSdata(params);
    else
        [Rtest,Stest] = DATAGENPP(Nbatches,params,datagenargs{:});
        [testData.R,testData.S] = longdata(Rtest,Stest);
        clear Rtest Stest
    end
    
    yvar = []; vvar = [];
    if usejava('desktop')
        setColors;
        figure(2014); clf; hold on;
        subplot(1,2,1); hold on;
        plotHandle(1) = plot(NaN,NaN);
        hold off;
        subplot(1,2,2); hold on;
        plotHandle(2) = plot(NaN,NaN);
        plotHandle(3) = plot(NaN,NaN);
        hold off;
    end
end



% NOW OVERWRITE WITH THE LOADED WEIGHTS
vishid = wts{i_rbm}(1:end-1,:);
hidbiases = wts{i_rbm}(end,:);
visbiases = wts{numRBMs*2-i_rbm+1}(end,:)';

% run 
EFH

% store the resulting weights!
wts{i_rbm} = [vishid; hidbiases];
wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases'];