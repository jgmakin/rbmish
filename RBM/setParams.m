function params = setParams(varargin)
% SETPARAMS     Sets parameters for sensory integration learning
%   SETPARAMS sets the (reused) parameters of DBN, DATAGEN, FK2link, etc.

%-------------------------------------------------------------------------%
% Revised: 09/29/16
%   -allowed for more general setting of parameters with variable input
%   arguments
% Revised: 08/25/16
%   -unified all the rEFHs that run on (linear) PPC data
%   -unified/functionized the learning-rate-update functions
% Revised: 12/09/13
%   -rationalized: put in case statements for all the different models
% Revised: 09/07/12
%   -fixed max and min calculations
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme and non-2D-stimuli.
% Revised: 7/14/10
%   -changed the grid size---while keeping the relative tuning acuity
%       constant.
% Revised: 7/7/10
%   -fixed posmax (should be measured from origin, not from left edge)
% Created: 7/5/10
%   by JGM
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%
% Consider restoring MODEL, but not as a field in params; rather, as a way
% of running your favorite combinations of different parameters.  Then this
% could be more conveniently called from the outside to, say, reproduce
% some of your old results.
%%%%%%%%%%%%%


% what computer?
[~,machine] = system('hostname');
params.machine = strtrim(machine);

% variable arguments
params.datatype     = defaulter('datatype','1Dinteg',varargin{:});
params.Ncdsteps     = defaulter('Ncdsteps',1,varargin{:});
params.Ncases       = defaulter('Ncases',40,varargin{:});         % 10
params.Nbatches     = defaulter('Nbatches',1000,varargin{:});     % ?
params.NepochsMax   = defaulter('NepochsMax',90,varargin{:});     % 50
SensoryUnitType     = defaulter('SensoryUnitType','Poisson',varargin{:});



% which model?
switch params.datatype
    
    case '2Dinteg'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % data generation
        % data generation
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 30;                      % number of neurons *** (2) ***
        params.NS = 'Joint-Angle';          % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = [-pi/2; pi/4];              % [-pi/4; 0];
        thmax = [pi/4; 3*pi/4];             % [pi/4; pi/2];
        armlengths = [12; 20];
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin];
        params.smax = [roboparams.posmax roboparams.thmax];
        params.roboparams = roboparams;
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis/2};
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}};
        %%%params.typeUnits = {'Bernoulli','BernoulliDropout'};
        %%%params.typeUnits = {'Binomial','Binomial'};
        
        % gains
        if any(strcmp(params.typeUnits{2},'BernoulliDropout'))
            params.gmin = [0.064, 0.064];
            params.gmax = [0.096, 0.096];
        else
            params.gmin = [12, 12];
            params.gmax = [18, 18];
        end
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % learning schedules
        params = setLearningSchedules(1000,250,250,'exp',params);
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
    
    case '2Dinteg_NonFlatPrior'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
      
        % data generation
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 30;                      % number of neurons *** (2) ***
        params.NS = 'Joint-Angle';          % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = [-pi/2; pi/4];              % [-pi/4; 0];
        thmax = [pi/4; 3*pi/4];             % [pi/4; pi/2];
        armlengths = [12; 20];
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin];
        params.smax = [roboparams.posmax roboparams.thmax];
        params.roboparams = roboparams;
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis/2};
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}};
        
        % gains
        params.gmin = [12, 12];
        params.gmax = [18, 18];
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % learning schedules
        params = setLearningSchedules(2000,500,500,'exp',params);
        params.NepochsMax = 10;
        %%%% change me
        p0.cov = 1.0e-03*[0.2467,0;0,0.1097];
        p0.mu = [-0.3927;1.5708];
        params.p0 = p0;
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params,'latentvarprior',p0));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
    case '2Dinteg_AllGains'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % data generation
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 30;                      % number of neurons *** (2) ***
        params.NS = 'Joint-Angle';          % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = [-pi/2; pi/4];              % [-pi/4; 0];
        thmax = [pi/4; 3*pi/4];             % [pi/4; pi/2];
        armlengths = [12; 20];
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin];
        params.smax = [roboparams.posmax roboparams.thmax];
        params.roboparams = roboparams;
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis/2};
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}};
        %%%params.typeUnits = {'Bernoulli','BernoulliDropout'};
        %%%params.typeUnits = {'Binomial','Binomial'};
        
        % gains
        if any(strcmp(params.typeUnits{2},'BernoulliDropout'))
            params.gmin = [0, 0];
            params.gmax = [0.1, 0.1];
        else
            params.gmin = [0, 0];
            params.gmax = [20, 20];
        end
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % learning schedules
        params = setLearningSchedules(250,50,50,'exp',params);
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(Ntestdata,...
            @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,setfield(setfield(params,'gmin',[2,2]),...
            'gmax',[18,18]))),params.getData,yrclass));
        
        
    case '2Dtwoarms'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        %%%% consider just calling both 'Joint-Angle'???
        params.mods{1} = 'Joint-Angle-Left';
        params.mods{2} = 'Joint-Angle-Right';
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 30;                      % number of neurons *** (2) ***
        params.NS = 'Joint-Angle-Left';    	% "neutral space"
        params.walls = 'clipping';
        params.gmin = [12 12];
        params.gmax = [18 18];
        
        % arm properties
        params.smin = [-pi/2 -pi/2; pi/4 pi/4];
        params.smax = [pi/4 pi/4; 3*pi/4 3*pi/4];
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis/2};
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}};
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % learning schedules
        params = setLearningSchedules(1000,250,250,'exp',params);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));        
        
        
        
    case '1Daddition'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.mods{3} = 'Gaze-Angle';
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.NS = 'Hand-Position';        % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = pi/6;
        thmax = 5*pi/6;
        armlengths = 12;
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin roboparams.eyemin];
        params.smax = [roboparams.posmax roboparams.thmax roboparams.eyemax];
        params.roboparams = roboparams;
        
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}};
        % params.typeUnits = {{'Bernoulli'},{'BernoulliDropout'}};
        Nmods = length(params.mods);
        if strcmp(params.typeUnits{2},'BernoulliDropout')
            params.gmin = [0.024 0.072 0.024];
            params.gmax = [0.036 0.108 0.036];
            params.N = 600;
            Nvis = Nmods*params.N^params.Ndims;
            params.numsUnits = {Nvis, 1200};
            params.NepochsMax = 25;
            params = setLearningSchedules(1,0.25,0.25,'exp',params);
            params.kw = 1e-7;                 % because of dropout!
            params.kvb = 1e-7;                %   but see choice for
            params.khb = 1e-7;                %   '1Dinteg' below
        else
            params.gmin = [4 12 4];
            params.gmax = [6 18 6];
            params.N = 60;                    % number of neurons *** (2) ***
            Nvis = Nmods*params.N^params.Ndims;
            params.numsUnits = {Nvis, 160};
            params = setLearningSchedules(250,60,60,'exp',params);
        end
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
    case '1Daddition_2Layer'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.mods{3} = 'Gaze-Angle';
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.NS = 'Hand-Position';        % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = pi/6;
        thmax = 5*pi/6;
        armlengths = 12;
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin roboparams.eyemin];
        params.smax = [roboparams.posmax roboparams.thmax roboparams.eyemax];
        params.roboparams = roboparams;
        
        params.typeUnits = {{'Poisson'}, {'Bernoulli'}, {'Bernoulli'}};
        Nmods = length(params.mods);
        params.gmin = [4 12 4];
        params.gmax = [6 18 6];
        params.N = 60;                    % number of neurons *** (2) ***
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, 160, 160};
        params = setLearningSchedules(250,60,60,'exp',params);
        
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
    case '1Dinteg'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.Ndims = 1;                   % encoded vars/neuron
        params.NS = 'Joint-Angle';
        params.walls = 'clipping';
        
        % arm properties
        thmin = pi/6;
        thmax = 5*pi/6;
        armlengths = 12;
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin];
        params.smax = [roboparams.posmax roboparams.thmax];
        params.roboparams = roboparams;
        
        
        % RBM units
        params.typeUnits = {{'Poisson'},{'Bernoulli'}};
        Nmods = length(params.mods);
        %%%params.typeUnits = {{'Bernoulli'},{'BernoulliDropout'}};
        if strcmp(params.typeUnits{2},'BernoulliDropout')
            params.gmin = [0.064, 0.064];
            params.gmax = [0.096, 0.096];
            params.N = 600;
            Nvis = Nmods*params.N^params.Ndims;
            params.numsUnits = {Nvis, Nvis/2};
            params.NepochsMax = 25;
            params = setLearningSchedules(1,0.25,0.25,'exp',params);
            %%% => k < 0.3125 for no osc
            params.kw = 1e-7;
            params.kvb = 1e-7;
            params.khb = 1e-7;
        else
            params.gmin = [12, 12];
            params.gmax = [18, 18];
            params.N = 60;
            Nvis = Nmods*params.N^params.Ndims;
            params.numsUnits = {Nvis, Nvis}; % Nvis/2};
            params.NepochsMax = 90;
            params = setLearningSchedules(500,120,120,'exp',params);
        end
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % set the data generation functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case '2Daddition'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.mods{3} = 'Gaze-Angle';
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 30;                      % number of neurons *** (2) ***
        params.gmin = [12 12 12];
        params.gmax = [18 18 18];
        params.NS = 'Joint-Angle';          % "neutral space"
        params.walls = 'clipping';
        
        % arm properties
        thmin = [-pi/2; pi/4];              % [-pi/4; 0];
        thmax = [pi/4; 3*pi/4];             % [pi/4; pi/2];
        armlengths = [12; 20];
        roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,...
            params.mods,params.NS);
        params.smin = [roboparams.posmin roboparams.thmin roboparams.eyemin];
        params.smax = [roboparams.posmax roboparams.thmax roboparams.eyemax];
        params.roboparams = roboparams;
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis};
        params.typeUnits = {{'Poisson'},{'Bernoulli'}};
        
        
        % learning schedules
        params = setLearningSchedules(500,120,120,'exp',params);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            params.NS,params));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'HierL2'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % load a saved model
        fileprefix = 'wts_2Dinteg_900';
        % oldvars = load(sprintf('%sRBMish/EFHs/%s_31-Dec-2016.mat',...
        %     getdir('data'),fileprefix),'wts','params');
        oldvars = load(sprintf('%sRBMish/EFHs/new/%s.mat',...
            getdir('data'),fileprefix),'wts','params');
        oldvars.params.smpls = 15; % to generate hidden-layer activities
        oldvars.params.machine = params.machine;
        
        % new (upper) EFH: *mostly* the same as the lower EFH
        reusedparams = {'Ts','Ndims','N','walls',...
            'C','respLength','margin','gridsize','granularity',...
            'Ncases','Nbatches','NepochsMax','Ncdsteps'};
        for fn = reusedparams
            params.(fn{1}) = oldvars.params.(fn{1});
        end
        
        % stimulus parameters (one modality)
        params.NS = 'Joint-Angle';      % say, the *left* arm
        params.mods = {'Joint-Angle'};
        modInds = strcmp(oldvars.params.mods,params.mods);
        params.smin = oldvars.params.smin(:,modInds);
        params.smax = oldvars.params.smax(:,modInds);
        params.gmin = oldvars.params.gmin(:,modInds);
        params.gmax = oldvars.params.gmax(:,modInds);
        
        % EFH parameters (two typesUnits{1})
        params.typeUnits = {{'Bernoulli','Poisson'},{'Bernoulli'}};
        Nhid0 = sum(oldvars.params.numsUnits{end});
        params.numsUnits = {[Nhid0, Nhid0], Nhid0};
        params.mw  = @(iii)([oldvars.params.mw(iii)/10,oldvars.params.mw(iii)]);
        params.mvb = @(iii)([oldvars.params.mvb(iii)/10,oldvars.params.mvb(iii)]);
        params.mhb = @(iii)(oldvars.params.mhb(iii));
        params.bw  = @(iii)([oldvars.params.bw(iii)/10,oldvars.params.bw(iii)]);
        params.bvb = @(iii)([oldvars.params.bvb(iii)/10,oldvars.params.bvb(iii)]);
        params.bhb = @(iii)(oldvars.params.bhb(iii));
        params.kw  = @(iii)([oldvars.params.kw(iii)/10,oldvars.params.kw(iii)]);
        params.kvb = @(iii)([oldvars.params.kvb(iii)/10,oldvars.params.kvb(iii)]);
        params.khb = @(iii)(oldvars.params.khb(iii));
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(...
            getLatentsMultisensoryIntegration(Nexamples,yrclass,...
            oldvars.params.NS,oldvars.params));
        params.getData = @(S,Q)(getDataHierL2(S,Q,oldvars.wts,...
            oldvars.params,params));
        params.testEFH = @(R,S,Q,wts)(testEFHPPC(R,S,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'LTI-PPC'
        
        % which learning algorithm?
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        %%% params.algorithm = defaulter('algorithm','TRBM',varargin{:});
        %%% params.algorithm = defaulter('algorithm','RTRBM',varargin{:});
        
        % include a spring force?
        SPRING = defaulter('SPRING',1,varargin{:});
        
        % which set of modalities?
        params.mods = defaulter('mods',{'Joint-Angle'},varargin{:});
        % params.mods = defaulter('mods',{'Joint-Angle','Angular-Velocity'},varargin{:});
        % params.mods = defaulter('mods',{'Joint-Angle','Efference-Copy'},varargin{:});
        % params.mods = defaulter('mods', {'Joint-Angle','Angular-Velocity','Efference-Copy'},varargin{:});
        
        % how many dimensions of positions?
        Ndims = defaulter('Ndims',1,varargin{:});
        
        % data parameters
        params.Ndims = Ndims;               % encoded vars/neuron
        params.walls = 'wrapping';
        
        % state properties
        zmax = [pi/3 0.8]; % second-order
        if any(strcmp(params.mods,'Efference-Copy')), zmax = [zmax 1.25]; end
        if ~SPRING, zmax(2) = 1.9; end
        zmax = repmat(zmax,[Ndims,1]); zmax = zmax(:);
        zmin = -zmax;
        params.dynamics = setDynamics(Ndims,params.mods,zmin(1),zmax(1));
        dt = params.dynamics.dt;
        Nmods = length(params.mods);
        params.smin = reshape(params.dynamics.C*zmin,[Ndims,Nmods]);
        params.smax = reshape(params.dynamics.C*zmax,[Ndims,Nmods]);
        
        % data generation
        switch SensoryUnitType              % gain, swin in gain
            case 'Poisson'
                params.N = 15;
                params.gmin = 6.4*ones(size(params.mods));
                params.gmax = 9.6*ones(size(params.mods));
                Nvis = Nmods*params.N^params.Ndims;
                switch Ndims
                    case 1, Nhid = 16*Nvis; % 19*Nvis;
                    case {2,3}, Nhid = 5*Nvis;
                end
            case 'Bernoulli'
                params.N = 15;
                params.gmin = 0.08*ones(size(params.mods));
                params.gmax = 0.08*ones(size(params.mods));
                Nwindows = 100;
                dt = dt/Nwindows; m = 5; b = 0.25; k = 3;
                params.dynamics.dt = dt;
                params.dynamics.A =...
                    [eye(Ndims)         dt*eye(Ndims);...
                    -k/m*dt*eye(Ndims)  -(b/m*dt-1)*eye(Ndims)];
                Nvis = Nmods*params.N^params.Ndims;
                Nhid = 4*Nvis;
            case 'StandardNormal'
                params.N = 1;
                sigmaSqY = 0.0012;
                params.gmin = 1/(sqrt(sigmaSqY))*ones(size(params.mods));
                params.gmax = 1/(sqrt(sigmaSqY))*ones(size(params.mods));
                params.dynamics.SigmaYX = sigmaSqY*eye(Nmods);
                params.walls = 'none'; %%% ??
                Nvis = Nmods*params.N^params.Ndims;
                switch Ndims
                    case 1, Nhid = 15*Nvis;
                    case {2,3}, Nhid = 5*Nvis;
                end
        end
        params.NS = 'Joint-Angle';
        
        % the "NoSpring" case
        % the "NoSpring" case
        if ~SPRING
            % to make the NoSpring case:
            k=0; b=0.25; m=5; dt=0.05;
            A = [eye(Ndims)             dt*eye(Ndims);...
                -k/m*dt*eye(Ndims)      -(b/m*dt-1)*eye(Ndims)];
            params.dynamics.A = A;
            params.dynamics.SigmaX = params.dynamics.SigmaX*50;
            params.dynamics.m = m;
            if strcmp(SensoryUnitType,'Bernoulli')
                Nhid = 1*Nvis;
            else
                Nhid = 15*Nvis;
            end
        end
        
        % RBM units, cont
        params.numsUnits = {Nvis, Nhid};
        params.typeUnits = {{SensoryUnitType},{'Bernoulli'}};
        
        % delayed?
        %%% params.delay = 2;
        
        % use different learning rates for Bernoulli-Bernoulli wts
        switch params.algorithm
            case {'TRBM','RTRBM'}
                % adjust learning parameters....
                params.Ncases = 100;
                params.Nbatches = 40000/params.Ncases;
                params.EACHBATCHISATRAJ = 1;
                params.Npretrain = 0; % s30;
                params.NepochsMax = 255;
                params.Ncdsteps = 25;
                params = setLearningSchedules(1500,1500,1500,'hyperbolic',params,150,150);
            case 'rEFH'
                params.EACHBATCHISATRAJ = 0;
                params.Npretrain = 0;
                if strcmp(SensoryUnitType,'Bernoulli')
                    % params = setLearningSchedules(50,12,12,'exp',params,1);
                    params = setLearningSchedules(1,0.25,0.25,'exp',params,1,0.25);
                    % %%% => k < 0.3125 for no osc
                    % params.kw = 1e-7;               % because of dropout!
                    % params.kvb = 1e-7;              %
                    % params.khb = 1e-7;
                    params.Ncases = 5;
                    params.Nbatches = 40000/params.Ncases;
                else
                    switch Ndims
                        case 1
                            if any(strcmp(params.mods,'Efference-Copy')),
                                params.NepochsMax = 255;
                                %params = setLearningSchedules(800,200,200,'logistic',params,50,12);
                                %params.NepochsMax = 90;
                                params.EACHBATCHISATRAJ = 1;
                                params.Ncases = 100;
                                params.Nbatches = 40000/params.Ncases;
                                params = setLearningSchedules(1500,1500,1500,'hyperbolic',params,150,150);
                            else
                                params.NepochsMax = 50;
                                params = setLearningSchedules(...
                                    800,200,200,'exp',params,50,12);
                            end
                        case {2,3}
                            params.NepochsMax = 1200;
                            params = setLearningSchedules(500,120,120,'logistic',params,50,12);
                    end
                end
            otherwise
                error('unknown training algorithm -- jgm');
        end
        
        % get the "grid parameters"
        [params.C,params.respLength,params.margin,params.gridsize,...
            params.granularity] = getGridParams(params.Ndims,params.N,...
            params.walls);
        
        % functions
        if params.EACHBATCHISATRAJ, T=params.Ncases; else T=params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(getLatentsLTI(...
            Nexamples,T,yrclass,params.NS,params,varargin{:}));
        params.getData = @(X,Q)(getDataPPC(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHPPC(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(40,1000,yrclass,params));
        
        
        
        
    case 'MCDexperiment'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods = {'Motion-Dots','ICMSpolar'};
        
        % just for this case
        dmin = 40;      % mm
        dmax = 115;     % mm
        Fmin = 5; % 100;     % Hz
        Fmax = 15; % 300;     % Hz
        Nelectrodes = 8;
        PDs = 2*pi*(0:Nelectrodes-1)/Nelectrodes - pi;
        
        params.experiment.Fmin = Fmin;
        params.experiment.Fmax = Fmax;
        params.experiment.dmin = dmin;
        params.experiment.dmax = dmax;
        params.experiment.Nelectrodes = Nelectrodes;
        params.experiment.PDs = PDs;
        
        electrodeFunc = @(th,d)(Fmax/dmax*(d.*(1+cos(th-PDs)))/2);
        
        params.experiment.electrodeFunc = electrodeFunc;
        params.experiment.neuronsPerElectrode = 30;
        
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 15;                      % 24^2 = 8*72
        params.gmin = [15 15];
        params.gmax = [15 15];
        params.NS = 'Motion-Dots';          % "neutral space"
        
        % arm properties
        cartmin = -[dmax; dmax];
        cartmax = +[dmax; dmax];
        polarmin = [0; 40];
        polarmax = [2*pi; 115];
        params.smin = [cartmin, polarmin];
        params.smax = [cartmax, polarmax];
        
        % RBM units
        Nvis = params.N^params.Ndims +...
            Nelectrodes*params.experiment.neuronsPerElectrode;
        params.numsUnits = {Nvis, floor(Nvis/2)};
        params.typeUnits = {{'Poisson'},{'Bernoulli'}};
        
        % learning schedules
        params = setLearningSchedules(500,120,120,'exp',params);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(getLatentsMCD(...
            Nexamples,yrclass,params));
        params.getData = @(S,Q)(getDataMCDexperiment(S,Q,params));
        params.testEFH = @(R,S,Q,wts)(testEFHPPC(R,S,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'MCDdarpa'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        % modality names
        params.mods = {'Motion-Dots','ICMS'};
        
        % distance min and max
        dmin = 40;      % mm
        dmax = 115;     % mm
        params.experiment.dmin = dmin;
        params.experiment.dmax = dmax;
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 15;                      % number of neurons *** (2) ***
        params.xpctGains = 15*ones(size(params.mods));
        params.swing = 1;                   % swing in gain (max=1)
        params.NS = 'Motion-Dots';
        
        % arm properties
        params.smin = [-dmax -dmax; -dmax -dmax];
        params.smax = [dmax dmax; dmax dmax];
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        params.numsUnits = {Nvis, Nvis/2};
        params.typeUnits = {{'Poisson'},{'Bernoulli'}};
        
        % learning schedules
        params = setLearningSchedules(500,120,120,'exp',params);
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(getLatentsMCD(...
            Nexamples,yrclass,params));
        params.getData = @(S,Q)(getDataPPC(S,Q,params));
        params.testEFH = @(R,S,Q,wts)(testEFHPPC(R,S,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'HHSreachData'
        %%%%%
        % FIX ME
        %%%%%
        
        
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        
        % modality names
        params.mods = {'Hand-Position','Hand-Velocity','Efference-Copy'};
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.N = 15;                      % number of neurons *** (2) ***
        % params.N = 19;
        params.gmin = [6.4 6.4 6.4];
        params.gmax = [9.6 9.6 9.6];
        params.NS = 'Hand-Position';
        
        % you got these by looking!
        params.xmin = [-55; -100];          % min(targets) - 10;
        params.xmax = [105; 60];            % max(targets) + 10;
        params.vmin = -[500; 500];
        params.vmax = [500; 500];
        params.ctrlmin = -[34700; 33200];
        params.ctrlmax =  [34700; 33200];
        params.smin = [params.xmin, params.vmin, params.ctrlmin];
        params.smax = [params.xmax, params.vmax, params.ctrlmax];
        
        % dynamical params (copied from getCenterOutReacher.m)
        Nstates = 4;
        Ninputs = 2;
        Noutputs = 4; %%% assume velocity is reported by the senses!
        params.dynamics.A =...
            [0.9405   -0.0317    0.0510    0.0009;...
            -0.0470    0.9443    0.0027    0.0533;...
            -1.8574   -0.9579    0.9365   -0.0099;...
            -1.5016   -1.8170    0.0726    1.0336];
        params.dynamics.muX = 1e-3*[0.7610; -0.0935; 25.9366; -3.4016];
        params.dynamics.SigmaX =...
            [0.0094    0.0032    0.3653    0.1277;...
            0.0032    0.0055    0.1248    0.2192;...
            0.3653    0.1248   14.6508    5.1692;...
            0.1277    0.2192    5.1692    8.9392];
        params.dynamics.muX0 = [21.5466; -18.0314];
        params.dynamics.muV0 = [-11.4596; 11.3564];
        % params.dynamics.SigmaX0 =...
        %  [0.3955    0.0639   -0.0922   -0.0848;...
        %  0.0639    0.3432   -0.0008   -0.0212;...
        %  -0.0922   -0.0008    0.7269    0.4551;...
        %  -0.0848   -0.0212    0.4551    0.4148];
        params.dynamics.SigmaX0 = [123.0287 -106.4835; -106.4835 96.7165];
        params.dynamics.SigmaV0 = 1e3*[6.2359 0.4711; 0.4711 8.5808];
        params.dynamics.dt = 0.05;
        params.dynamics.C = eye(Noutputs,Nstates);
        params.dynamics.G = [0, 0; 0, 0; 0.01, 0; 0 0.01];
        params.dynamics.H = eye(Ninputs);
        %%%%%
        % for certain purposes, you still need the parameters that control
        % the "evolution" of the control.  But you may be able to get away
        % with skipping them if you just fit the LDS with EM or fullyObs,
        % rather than using "the true params."
        %%%%%
        params.walls = 'none'; %%% needs to be set, but doesn't matter what
        
        % RBM units
        Nmods = length(params.mods);
        Nvis = Nmods*params.N^params.Ndims;
        Nhid = 320;
        params.numsUnits = {Nvis, Nhid};
        params.typeUnits = {{'Poisson'},{'Bernoulli'}};
        
        % learning parameters
        params.Ncases = 200;
        params.Nbatches = 40000/params.Ncases;
        params.EACHBATCHISATRAJ = 1;
        params.Npretrain = 30;
        params.NepochsMax = 255;
        params.Ncdsteps = 25;
        params = setLearningSchedules(1500,750,750,'hyperbolic',params,150,75);
        
        % params.NepochsMax = 25000;
        % logisticGrowth = @(x0,ii)(1000/(1 + exp(-ii/8 + 7.5)) + 0.5)*x0*Ts^2;
        % mw0 = 500;
        % params.mw = @(iEp)(logisticGrowth(mw0,iEp));
        % params.mvb = @(iEp)(logisticGrowth(120,iEp));
        % params.mhb = @(iEp)(logisticGrowth(120,iEp));
        % fprintf('Using the SIGMOIDAL learning-rate adjustment scheme!!\n');
        % params.b = @(iEp)(0.5*mw0*Ts);
        
        if params.EACHBATCHISATRAJ, T=params.Ncases; else T=params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(getLatentsHHS(...
            Nexamples,T,yrclass,params.NS,params,varargin{:}));
        params.getData = @(S,Q)(getDataPPC(S,Q,params));
        params.testEFH = @(R,S,Q,wts)(testEFHPPC(R,S,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(40,1000,yrclass,params));
        
        
        
        
    case 'bouncingballs'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        % params.algorithm = defaulter('algorithm','TRBM',varargin{:});
        % params.algorithm = defaulter('algorithm','RTRBM',varargin{:});
        % params.algorithm = defaulter('algorithm','EFH',varargin{:});
        PAPER = '[Sutskever2013]';
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        
        switch PAPER
            case '[Sutskever2013]'
                params.N = 30;
                Nhid = 400;
                params.Ncases = 100;
                params.Nbatches = 400;
                params.NepochsMax = 255;
            case '[Boulanger-Lewandowski2012]'
                params.N = 15;
                Nhid = 300;
                params.Ncases = 100;
                params.Nbatches = 400;
                params.NepochsMax = 255;
            case '[Mittelman2014]'
                params.N = 30;
                Nhid = 2500;
                params.Ncases = 100;
                params.Nbatches = 100;
                params.NepochsMax = 1000;
            otherwise
                params.N = 30;
                Nhid = 625;
                params.Ncases = 100;
                params.Nbatches = 400;
                params.NepochsMax = 255;
        end
        
        % arm properties
        params.smin = [0;0];
        params.smax = [params.N;params.N];
        
        % modality names
        params.mods = {'balls'};
        params.walls = 'bouncing';
        %%% not actually used
        
        % RBM units
        Nsensory = params.N^params.Ndims;
        params.numsUnits = {Nsensory, Nhid};
        params.typeUnits = {{'Bernoulli'},{'Bernoulli'}};
        
        % learning parameters
        params.EACHBATCHISATRAJ = 1;
        params.Npretrain = 30;
        params.Ncdsteps = 25;
        params = setLearningSchedules(50,50,50,'hyperbolic',params,50,50);
        %params = setLearningSchedules(25,25,25,'hyperbolic',params,25,25);
        
        % balls
        params.balls.N = 3;                             % IS: 3
        params.balls.boundingbox = 10;                  % IS: 10
        params.balls.r = 1.2*ones(1,params.balls.N);    % IS: 1.2
        params.balls.speed = 0.5;                       % IS: 0.5
        
        if params.EACHBATCHISATRAJ, T=params.Ncases; else T=params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(getLatentsBouncingBalls(...
            Nexamples,yrclass,T,params.balls,ones(1,params.balls.N),varargin{:}));
        params.getData = @(X,Q)(getDataBouncingBalls(X,Q,params.balls,params.N));
        params.testEFH = @(R,X,Q,wts)(testEFHNextFrameError(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(400,100,yrclass,params));
        
        
        
        
    case 'ErlangMixtureToy'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        
        
        % RBM units
        Nsamples = 1; %180; % can be done with 8
        Ncats = 4;
        params.numsUnits = {Nsamples*2, Ncats-1};
        params.typeUnits = {{'Erlang'},{'Categorical'}};
        %%% Bernoulli also works fine
        
        % Erlang parameters
        %         params.shapeparams = [1;16;7];          %   1   16	7
        %         params.scaleparams = [0.05;0.01;0.03];  %   20  100 33
        %         ps = [10;3;3];
        %         params.mixingproportions = ps(1:end-1)/sum(ps); % [p1,...,p_{n-1}]
        %
        %         params.shapeparams = [2; 10];
        %         params.scaleparams = [0.04; 0.035];
        %         params.mixingproportions = 0.8;
        %
        %         params.shapeparams = [2.5; 5.5];
        %         params.scaleparams = [0.027; 0.058];
        %         params.mixingproportions = 0.73;
        params.shapeparams = [2.5; 5.5];
        params.scaleparams = [0.027; 0.054];
        params.mixingproportions = 0.77;
        
        % learning (hyper-)parameters
        params.NepochsMax = 25;
        params.Ncases = 50;
        params.Nbatches = 40000/params.Ncases;
        %%%% surely the varargin here shouldn't be the default, 1/10!??!
        params = setLearningSchedules(8,8,10,'hyperbolic',params);
        %params = setLearningSchedules(15,15,25,'hyperbolic',params);
        %params = setLearningSchedules(60,60,100,'hyperbolic',params);
        
        %%% seems like: as the ratio of visibles to hiddens *increases*,
        %%% you have to turn down the momentum accordingly....
        
        % functions
        params.getLatents = @(Nexamples,yrclass)(getLatentsClasses(...
            Nexamples,yrclass,params));
        params.getData = @(X,Q)(getDataErlangMixture(X,Q,params));
        params.testEFH = @(R,X,Q,wts)(testEFHErlangMixture(R,X,Q,wts,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'ECcoherences'
        params.algorithm = defaulter('algorithm','EFH',varargin{:});
        % params.algorithm = defaulter('algorithm','rEFH',varargin{:});
                
        % RBM units
        Nsamples = 1; % 80; %%% 180*12/60 = 36 minutes
        Ncats = 4;
        %%%
        
        % units
        params.numsUnits = {Nsamples*2, Ncats-1};
        params.typeUnits = {{'Erlang'},{'Categorical'}};
        %params.numsUnits = {Nsamples*2, log2(Ncats)};
        %params.typeUnits = {{'Erlang'},{'Categorical'}};
        %params.numsUnits = {Nsamples, log2(Ncats)};
        %params.typeUnits = {{'GammaFixedScale'},{'Categorical'}};
        %params.numsUnits = {Nsamples, log2(Ncats)};
        %params.typeUnits = {'GammaFixedScale'},{'Bernoulli'}};
        %params.numsUnits = {Nsamples*2, log2(Ncats)};
        %params.typeUnits = {{'Erlang'},{'Bernoulli'}};
        
        % Erlang parameters
        %params.shapeparams = [0;0];         % we don't know these!
        %params.scaleparams = [0;0];         %
        params.scaleparams = 0.05;
        params.mixingproportions = 0.5;     %
        
        % data params
        params.Nsperwinstep = 1;
        params.Nsperwindow = 1;
        params.subj         = defaulter('subject','EC108',varargin{:});
        params.trainsuffix  = defaulter('trainsuffix','_8337e405',varargin{:});
        
        % learning (hyper-)parameters
        if strcmp(params.algorithm,'EFH')
            %%%% these learning rates are just made up....
            params.NepochsMax = 40;
            params.Ncases = 50;
            params.Nbatches = 40000/params.Ncases;
            params = setLearningSchedules(30,30,60,'hyperbolic',params,30,30);
            params.getLatents = @(Nexamples,yrclass,varargin)(...
                getLatentsCoherences(Nexamples,yrclass,params));
        else
            params.NepochsMax = 40;
            params.Ncases = 50;
            params.Nbatches = 40000/params.Ncases;
            params.EACHBATCHISATRAJ = 0;
            params.Npretrain = 0;
            if params.EACHBATCHISATRAJ
                T=params.Ncases;
            else
                T=params.Nbatches;
            end
            params = setLearningSchedules(30,30,60,'hyperbolic',params,30,30);
            params.getLatents = @(Nexamples,yrclass,varargin)(...
                getLatentsCoherences(Nexamples,yrclass,params,...
                'sequencelength',T,varargin{:}));
        end
        params.getData = @(S,Q)(getDataCoherences(S,Q,params.typeUnits{1}(1)));
        params.testEFH = @(R,X,Q,wts)(testEFHcoherences(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(160,250,yrclass,params));
        Ntestdata = 40000;
        params.getTestData = @(yrclass)(generateData(...
            Ntestdata,params.getLatents,params.getData,yrclass));
        
        
        
        
    case 'ECoG'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.N = 20; %124;                     % number of neurons *** (2) ***
        
        % arm properties
        params.smin = 0;
        params.smax = params.N;
        
        
        % RBM units
        Nsensory = params.N^params.Ndims;
        Nhid = 400;
        params.numsUnits = {Nsensory, Nhid};
        params.typeUnits = {{'StandardNormal'},{'Bernoulli'}};
        
        % use different learning rates for Bernoulli-Bernoulli wts
        params.EACHBATCHISATRAJ = 0;
        params.Npretrain = 0;
        params.NepochsMax = 40;
        params.Ncdsteps = 25;
        params.Ncases = 20;
        params.Nbatches = 40000/params.Ncases;
        params = setLearningSchedules(800,200,200,'exp',params,130,30);
        %%% should the *hidden biases* (hb) be updated fast?
        
        % functions
        if params.EACHBATCHISATRAJ, T=params.Ncases; else T=params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(...
            getLatentsNull(Nexamples,yrclass,T,varargin{:}));
        params.getData = @(S,Q)(getDataECoG(S,Q,params.machine));
        params.testEFH = @(R,X,Q,wts)(testEFHNextFrameError(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(40,1000,yrclass,params));
        
        
        
        
    case 'spikecounts'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        datafile = defaulter('datafile',...
            'Indy_datafiles/spikes_and_kinematics_20160407_02.mat',...
            varargin{:});
        params.datafile = datafile;
        params.mods = defaulter('mods',{'M1S1'},varargin{:});
        params.Nmsperbin = defaulter('Nmsperbin',64,varargin{:});
        params.trainingtime = defaulter('trainingtime',320,varargin{:});
        % 320 seconds; also can choose 'half'.
        params.fraction = defaulter('fraction',1,varargin{:});
        
        if any(strcmp(datafile(1:4),{'Indy','Loco','Jack'}))
            params.Fs = 250;
        elseif strcmp(datafile(1:6),'Chewie')
            params.Fs = 1000;
        else
            error('unexpected BMI data file -- jgm');
        end

        
        % dynamical params
        params.typeUnits = {{SensoryUnitType},{'Bernoulli'}};
        
        % RBM units
        [~,Q] = getLatentsBMI([],'double',[],'train',params,...
            'sequencelength','singlesequence');
        [NsamplesTrain,Nsensory] = size(Q.R);
        Nhid = Nsensory*4;
        params.numsUnits = {Nsensory, Nhid};
        
        % batch division
        params.EACHBATCHISATRAJ = any(strcmp(params.algorithm,{'RTRBM','TRBM'}));
        params.Npretrain = 0;
        params.NepochsMax = max(floor(NsamplesTrain/250),20);

        % make sure trajs are long enough 
        Ntrainingsamples = params.trainingtime/(params.Nmsperbin/1000);
        if params.EACHBATCHISATRAJ
            params.Ncases = min(40,Ntrainingsamples);
            params.Nbatches = floor(40000/params.Ncases);
        else
            params.Nbatches = min(1000,Ntrainingsamples);
            params.Ncases = floor(40000/params.Nbatches);
        end
        
        % sparsity? (0.20,0.95,0.05)
        sparse.cost = 0.20;
        sparse.phidNewFrac = 0.95;
        sparse.phidTarget = 0.05;
        params.sparse = sparse;
        
        % learning rates
        muR = mean(Q.R(:));
        muZ = 0.09;             % \approx, because phidTarget = 0.05;
        if strcmp(params.typeUnits{1},'StandardNormal')
            mw0 = 200;
        else
            mw0 = mean([muR*Nsensory,muZ*Nhid])*10000/NsamplesTrain;
        end
        mr0 = 110;
        wt2bias = 0.55;         %
        params = setLearningSchedules(mw0,mw0*wt2bias,mw0*wt2bias,'exp',...
            params,mr0,mr0*wt2bias);
        
        % functions
        if params.EACHBATCHISATRAJ, T = params.Ncases; else T = params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(getLatentsBMI(...
            Nexamples,yrclass,T,'train',params,varargin{:}));
        params.getData = @(X,Q)(deal(Q.R, rmfield(Q,'R')));
        params.testEFH = @(R,X,Q,wts)(testEFHBMI(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(...
            1,'singlesequence',yrclass,setfield(params,'getLatents',...
            @(Nexamples,yrclass,varargin)(getLatentsBMI(...
            [],yrclass,[],'test',params,varargin{:})))));
        
        
    case 'MoCap'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        
        % training parameter        
        params.EACHBATCHISATRAJ = 0;
        params.Npretrain = 0;
        Nsensory = 49; % params.N^params.Ndims;
        
        %{
        Nhid = 80;
        params.NepochsMax = 80;
        params.Nbatches = 100;                  % => sequence length=100
        params.Ncases = 40000/params.Nbatches;
        params.numsUnits = {Nsensory, Nhid};
        params = setLearningSchedules(100,50,50,'exp',params,16,8);
        %}
        
        %%{
            Nhid = 200;
            params.NepochsMax = 150;
            params.Nbatches = 100;                  % => sequence length=100
            params.Ncases = 40000/params.Nbatches;
            params.numsUnits = {Nsensory, Nhid};
            params = setLearningSchedules(100,50,50,'exp',params,16,8);
        %%}
        
        params.typeUnits = {{'StandardNormal'},{'Bernoulli'}};
        path(path,genpath('C:\#code\MHMUBLV'));
        
        % functions
        if params.EACHBATCHISATRAJ, T = params.Ncases; else T = params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(...
            getLatentsNull(Nexamples,yrclass,T,varargin{:}));
        params.getData = @(S,Q)(getDataMoCap(S,Q,'train'));
        params.testEFH = @(R,X,Q,wts)(testEFHNextFrameError(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(800,50,yrclass,...
            setfield(params,'getData',@(S,Q)(getDataMoCap(S,Q,'test')))));
        
        
        
    case 'polyphonicmusic'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        
        % modality names
        params.mods = {'none'};
        
        
        % RBM units
        Nsensory = 88;
        Nhid = 120;  %%% first guess
        % Nhid = 400;
        params.numsUnits = {Nsensory, Nhid};
        params.typeUnits = {{'Bernoulli'},{'Bernoulli'}};
        
        % use different learning rates for Bernoulli-Bernoulli wts
        params.Ncases = 100;
        params.Nbatches = 400;
        params.EACHBATCHISATRAJ = 0;
        params.Npretrain = 0;
        params.NepochsMax = 80;
        params = setLearningSchedules(100,50,50,'exp',params,100,50);
        % params = setLearningSchedules(800,200,200,'exp',params,800,200);
        
        % which set of music
        % params.database = 'Nottingham';
        params.database = 'JSB Chorales';
        % params.database = 'MuseData';
        % params.database = 'Piano-midi';
        
        % functions
        if params.EACHBATCHISATRAJ, T = params.Ncases; else T = params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(...
            getLatentsNull(Nexamples,yrclass,T,varargin{:}));
        params.getData = @(S,Q)(getDataPolyphonicMusic(S,Q,...
            params.database,'train'));
        params.testEFH = @(R,X,Q,wts)(testEFHNextFrameError(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(800,50,yrclass,...
            setfield(params,'getData',@(S,Q)(getDataPolyphonicMusic(S,Q,...
            params.database,'test')))));
        
        
        
        
    case 'HDFRE'
        params.algorithm = defaulter('algorithm','rEFH',varargin{:});
        
        SensoryUnitType = 'Categorical';
        % SensoryUnitType = 'Bernoulli';
        Nchars = 81;
        switch SensoryUnitType
            case 'Categorical'
                Nsensory = Nchars-1;
            case 'Bernoulli'
                Nsensory = ceil(log2(Nchars));
            otherwise
                error('unexpected type of units for %s -- jgm',params.datatype)
        end
        Nhid = 100;
        params.numsUnits = {Nsensory,Nhid};
        params.typeUnits = {{SensoryUnitType},{'Bernoulli'}};
        params.mods = {'none'};
        
        params.EACHBATCHISATRAJ = 1;
        params.Ncases = 400;
        params.Nbatches = 100;
        params.Npretrain = 30;
        params.NepochsMax = 25;
        %params = setLearningSchedules(100,50,50,'exp',params,100,50);
        %params = setLearningSchedules(10,10,10,'exp',params,100,50);
        params = setLearningSchedules(100,100,100,'hyperbolic',params,100,100);
        
        % character set
        gibbonsymbols = ' !"&''()*,-./0123456789:;=?[]_';
        lowercaseletters = 'abcdefghijklmnopqrstuvwxyz';
        capitalletters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        params.charset = [gibbonsymbols,capitalletters,lowercaseletters];
        
        % functions
        if params.EACHBATCHISATRAJ, T = params.Ncases; else T = params.Nbatches; end
        params.getLatents = @(Nexamples,yrclass,varargin)(...
            getLatentsNull(Nexamples,yrclass,T,varargin{:}));
        params.getData = @(S,Q)(getDataHDFRE(S,Q,params.charset,params.typeUnits{1}{1}));
        params.testEFH = @(R,X,Q,wts)(testEFHNextFrameError(R,X,Q,wts,params));
        params.getTestData = @(yrclass)(dynamicalDataWrapper(40,1000,yrclass,params));
        
        
        
        
    otherwise
        
        error('unrecognized model -- jgm');
        
end


% universal PPC params


% setColors
params.VIScolor = [1 0 1];                  % magenta
params.PROPcolor = [1 0.5 0];               % orange
params.EYEcolor = [0 0 1];                  % blue
params.EFHcolor = [0 1 0];                  % green
params.OPTcolor = [0 0 0];                  % black


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [C,respLength,margin,gridsize,granularity] = getGridParams(Ndims,N,walls)

FWHM = 1/6;                             % deg (Deneve01: "45-75deg")
c = FWHM/(2*sqrt(2*log(2)));            % std. dev. in degrees
C = eye(Ndims)*c^2;
%%% essentially, one std of the tuning curve (in one direction only) is
%%% about 7% of the total grid; so the 1-std (in both directions) coverage
%%% of each tuning curve is 14% of the total area.

% grid parameters
respLength = 1;                         % normalized to 1  *** (3) ***
margin = 4*c;                           % 99.99%           *** (4) ***
if strcmp(walls,'wrapping'), margin = 0; end
gridsize = margin + respLength + margin;
granularity = N/gridsize;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [S,Q] = getLatentsNull(Nexamples,yrclass,T,varargin)

S           = zeros(Nexamples,0,1,yrclass);
Q.T         = defaulter('sequencelength',T,varargin{:});
Q.G         = [];
Q.states    = [];
Q.restarts  = [];
%%% may actually want to populate restarts, based on Ntraj....

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [R,Q] = getDataNull(X,Q,visDstrbs)
% USAGE:
%   R = getDataNull(S,Q,params.typeUnits{1}{end});

switch visDstrbs
    case 'StandardNormal'
        R = (Q.R - mean(Q.R))/chol(cov(Q.R));
        %%%%
        % consider just z-scoring, after the fashion of Suts
        %%%%
    case 'Bernoulli'
        R = Q.R;
        R(R>1) = 1;
    otherwise % probably 'Poisson'
        R = Q.R;
end

% strip R out of Q (no longer needed)
Q = rmfield(Q,'R');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [Rtest,Xtest,Qtest] = dynamicalDataWrapper(Ntraj,T,dataclass,params)


if isfield(params,'mods')&&isfield(params,'Ndims')
    if (params.Ndims == 1)&&~isfield(params,'delay')&&length(params.mods)==1
        load([getdir('data'),'RBMish/testdata/data_1D_LTI-PPC.mat']);
        if checkGPUavailability
            Rtest = gpuArray(Rtest);
            Xtest = gpuArray(Xtest);
        end
    elseif (params.Ndims==1)&&~isfield(params,'delay')&&...
            (length(params.mods)==2)&&any(strcmp(params.mods,'Efference-Copy'))
        load([getdir('data'),'RBMish/testdata/data_1D_LTI-PPC_withEC.mat']);
        if checkGPUavailability
            Rtest = gpuArray(Rtest);
            Xtest = gpuArray(Xtest);
        end
    else
        getLatents  = @(Nexamples,yrclass)(params.getLatents(Nexamples,...
            yrclass,'sequencelength',T));
        getData     = params.getData;
        [Rtest,Xtest,Qtest] = generateData(Ntraj*T,getLatents,getData,dataclass);
    end
else
    getLatents  = @(Nexamples,yrclass)(params.getLatents(Nexamples,...
        yrclass,'sequencelength',T));
    getData     = params.getData;
    [Rtest,Xtest,Qtest] = generateData(Ntraj*T,getLatents,getData,dataclass);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [X,Q] = getLatentsClasses(Nexamples,yrclass,params)

ps = params.mixingproportions;
if strcmp(yrclass,'gpuArray'), ps = gpuArray(ps); end
X = sampleT(repmat(ps',[Nexamples,1]),{'Categorical'},length(ps),params);
Q = [];


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function C = getNewObservationMatrix(H,SigmaYX)
% Confusing: You want to construct the data such that conditional variance
% of R, the neural response, is I; but the conditional variance of Y, what
% (e.g.) the KF sees, is SigmaYX---all while keeping the underlying state Z
% unchanged.  Recall also that the stimulus S = C*Z.  Hence, suppose
%   R ~ N(C*Z,I)
%   Y = M*R
% for some undetermined C and M.  How should we choose these so
%   E[Y|Z] = H*Z
%   Cov[Y|Z] = SigmaYX
% for some (determined) SigmaYX and H?  Well,
%   E[Y|Z] = M*C*Z = M*S
%   Cov[Y|Z] = M*Cov[R|Z]*M' = M*M'.
% Hence,
%   M := chol(SigmaYX)'
%   C := inv(M)*H,
% so indeed we have:
%   R ~ N(S,I)
%   Y ~ N(H*Z,SigmaYX)

M = chol(SigmaYX)';
C = M\H;

end
%-------------------------------------------------------------------------%


% *** (1) ***
% the reachable workspace is hard to state in general, so I've hard-coded
% in the relevant params---see drawing in your notes.  Ought to be changed

% *** (2) ***
% cf. Ma et al.: "p(s|r) converges to a Gaussian as the number of neurons
%  increases."

% *** (3) ***
% Deneve01 has her neurons cover 45-75 degrees---so use 60 degrees.  Now,
% I think the space she's working with spans 360 degrees, so the FWHM for
% each of her neurons is 1/6 (60/360) the total space.  You want to keep
% the same *relative* precision, so your neurons will have FWHM, not 60,
% but 1/6 the total space---which you have normalized to unity.
% ---But in fact [06/06/11], in the 2D case this doesn't even matter! since
% the tuning-curve widths are irrelevant there.

% *** (4) ***
% The margin should include at least 99.99% of the Gaussian bubbles; so for
%  bubbles at the edge of the respLength, there needs to be a 4 stddev
%  margin.

% *** (5) ***
% The control limits were were determined empirically: 0.06 is
%   approximately 3 stddevs out from the mean control (0) that was used to
%   generate Lissajous curves *with no limit on the magnitude of the
%   control*.  The limits were subsequently put into
%   generateLissajousControls.m.

% *** (6) ***
% Rationale for the IC priors:
% Total number of trajectories (and therefore ICs):
%   NepochsMax*Ncases = 90*40 = 3600
% Want E[number of rejections] < 1
% From lab notes:
%   (1) v0_rare = sqrt(umax/dt/A*(A^2 - (x0_rare - b)^2)).
% If x0_rare  = x0_{>2SD} and v0_rare = v0_{>3SD} then
%   (2) E[#rej] = (1-0.997)*(1-0.954)*3600 = 0.5 < 1 (ok)
% And x0 needs to be kept far away from the edges--say 5 SDs.
% Then:
%   (3) x0minusb_rare = 4*(xmax - xmin)/10; % (2SD each way)
%   (4) v0_rare = sqrt(0.06/dt./A.*(A.^2 - x0minusb_rare.^2))
% And 1 SD of v0 is:
%   (5) v0_{1SD} = v0_rare/3 = [0.5447, 0.4448]
% See assembleLissajousTensor for more details






%-------------------------------------------------------------------------%
% code cleanup, July 2014:
% (1) Check/rehabilitate the following (not a complete list):
%       standardCond2.m
%       gainmods.m
%       shiftTest.m
%       [yr adaptation files]
%
% (2) Try to eliminate the following:
%       ntrlJacobian.m
%       estError.m
%       estGather.m
%       estStatsCorePP.m
%       decode.m
%       decoder.m
%       marginalErrors.m
%       marginalErrorsPP.m
%       PPCinfo.m
%       dispErrCovs.m
%       covInCalc.m
%       PPCinputStats.m
%       SICE.m
%       tuner.m [appears not to be in use]
%
% (3) See if the following are still useful:
%       testStim.m
%       condInf.m
%       condInfPP.m
%       condInfCondErrs.m
%       confabulate.m
%       getNoiseFloorPP.m
%       optMeanCheck.m
%       rbmPP.m
%
% (N) back up all the "final wts" to zamfir and mushroom
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Good settings for learning "mass":
%
% Poisson-Bernoulli wts: 500
% Bernoulli biases: 120
% Poisson biases: 120
%
% Poisson-Bernoulli wts, hyperbolic learning: 1500
%
% RTRBM, TRBM, rEFH (Poisson-Bernoulli wts): 120
%
% Bernoulli-Bernoulli wts: 100
%
% Erlang-Categorical wts: 15 (but maybe b/c the weight matrix is so small)
% Categorical biases: 25
%-------------------------------------------------------------------------%
