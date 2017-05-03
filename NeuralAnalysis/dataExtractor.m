function TuningData = dataExtractor(DATE, Epoch, Align, DOPLOT)

DATADIR = 'C:\#DATA\Dmitri_datafiles';

% init
if(~exist('Align','var')),  Align='Variable';      end
if(~exist('Epoch','var')),  Epoch={'EyeOn' 'Go'};  end
if(~exist('DOPLOT','var')), DOPLOT=1;  end

% General Parameters
TrgType = 'TrgAng';  % 'TrgAng' or 'MvVecAng'


%% Load data
load([DATADIR,'/',DATE,'_behav.mat'],...
    'TuningTrialIndices_by_Target','ExperimentTrialIndices_by_Target',...
    'TuningTrialIndices','ExperimentTrialIndices');
load([DATADIR,'/',DATE,'_bparams.mat'],'All_TrgAngs','All_InitAngs',...
    'All_MvVectAngs','All_TrgAngs');
load([DATADIR,'/',DATE,'_neural.mat'],'NEURAL','ALL_TUNDATA','EPOCH_RATES');
% eval(sprintf('load %s/%s_mvtundata.mat MV_TUNDATA;', DATADIR,DATE));


%% Get EPOCH_RATES of interest
% find the element of EPOCH_RATES that matches epoch & alignment
epoch_starts = cellstr(strvcat(EPOCH_RATES(:).epoch_start));
epoch_ends = cellstr(strvcat(EPOCH_RATES(:).epoch_end));
align_list = cellstr(strvcat(EPOCH_RATES(:).align));
start_idx = strcmp(Epoch{1},epoch_starts);
end_idx = strcmp(Epoch{2},epoch_ends);
align_idx = strcmp(Align,align_list);
rates_idx = find(start_idx & end_idx & align_idx);
All_AllRates = EPOCH_RATES(rates_idx).rates;

%% Collect Data
tunTRG = All_TrgAngs(TuningTrialIndices);
tunIA  = All_InitAngs(TuningTrialIndices);
tunIA  = mod(tunIA-tunTRG+180,360)+tunTRG-180;
tunRATES = All_AllRates(TuningTrialIndices,:);

TuningData.tunTRG   = tunTRG;
TuningData.tunIA    = tunIA;
TuningData.tunRATES = tunRATES;

expTRG = All_TrgAngs(ExperimentTrialIndices);
expIA  = All_InitAngs(ExperimentTrialIndices);
expIA  = mod(expIA-expTRG+180,360)+expTRG-180;
expRATES = All_AllRates(ExperimentTrialIndices,:);

TuningData.expTRG   = expTRG;
TuningData.expIA    = expIA;
TuningData.expRATES = expRATES;

tunN = length(tunIA);
expN = length(expIA);


end