function varargout = DATAGENPP(Nbatches,params,varargin)
% DATAGENPP    Generate input vectors for the EFH/DBN
%
%   [R,S] = DATAGENPP(Nbatches,params,varargin)
%   [R,S,G] = DATAGENPP(Nbatches,params,varargin)
%   [R,S,G,Q] = DATAGENPP(Nbatches,params,varargin)
%
%   DATAGENPP randomly generates stimuli (S), their "neural responses" (R)
%   for a population of Gaussian-tuned neurons that smoothly tile the
%   space, and (if requested) the population gains (G).  The default dstrb 
%   of stimuli is uniform in the space of joint angles.
%
%   For other options, see the function specialvars below.
%
%   NB the sizes of the outputs:
%
%       R: (Nexamples x Nmods*N^Ndims)
%       S: (Nexamples x Ndims x Nmods)
%
%
%   VARIABLE ARGUMENTS:
%        string |   var   | range  | fxn                     | default val.
%   -----------------------------------------------------------------------
%      'swing'  | params. | >= 0   | overall max. fractional | params.swing
%               |  swing  |        |    gain change +/-;     |
%  'correlation'| prob    | [0,1]  | fraction of trials for  | 1
%               |         |        |    which x = f(th)      |
%  'vis/prop/eye| biases  | reals  | intermodal discrepancy  | zeros
%       bias'   |         |        |                         |
%  'vis/prop/eye| gains   | >=0    | mean max. spike counts  | params.g
%       gain'   |         |        |                         |
%     'prior'   |.mu,.cov | reals  | forces Gaussian prior   | uniform/
%               |         |        |    over the stimulus    | params.prior
%     'source'  | SOURCE  |vis,prop| which var gets sampled  | 'prop'
%     'verbose' | VERBOSE | {0,1}  | fprintf stuff?          | 1
%   'prevparams'|paramsOLD| params | for hierarchical nets   | void
%  'prevweights'| wtsOLD  |  wts   | for hierarchical nets   | void
%   'dbndepth'  | iRBM    |  Z+    | current layer of DBN    | void
%    'dbnwts'   |  wts    |  wts   | produce training data   | void
%               |         |        |  for deep layers of DBN |
%   'stimuli'   |  stims  |   S    | use existing stimuli    | void
%     'gains'   |  gains  |   G    | use existing gains      | void
% 'deadneurons' | deadInds| vector | zero out responses of   | void
%               |         | of Z+  |  these units            |


%-------------------------------------------------------------------------%
% Revised: 01/06/15 (JGM)
%   -added 'gains' as a varargin, and then used this hook in the toroidal
%   encoding (eliminating the need to set params.swing to 0 and then back
%   to its old value).
%   -updated the help section for the variable arguments
% Revised: 07/02/14 (JGM)
%   -re-wrote references to FK2link and IK2link to make use of their
%   vectorized forms
%   -this eliminated the last parfor loop!
% Revised: 06/13/14 (JGM)
%   -rewrote encodeStimuli in terms of case statements for different models
%   -added cases for 'MCDexperiment' throughout
% Revised: 05/06/14 (JGM)
%   -changed encodeStims to use vectorized respfxn.m, rather than a parfor
%   loop---so the only parallel loop left is in getStimuliCore.m.
% Revised: 12/18/13 (JGM)
%   -added case for controlled trajectories (tracking of Lissajous curves)
% Revised: 12/17/13 (JGM)
%   -added case for hierarchical data
%   -added internal data structure Q
% Revised: 12/10/13 (JGM)
%   -radically re-worked, especially getStimuli and encodeStimuli (formerly
%   getsourcevars and encodestim), in preparation for the addition of data
%   generated according to a controlled dynamical system
%   -output variable S (formerly X) now has dimensions Ncases x Ndims x 
%   Nmods x Nbatches, rather than Ncases x Ndims*Nmods x Nbatches.  This
%   matches the change in estStatsCorePP.m.
% Revised: 07/01/13 (JGM)
%   -changed getsourcevars to return a structure "datainfo" for the case of
%   a dynamic stimulus.
% Revised: 06/27/13 -BKD 
%   -for dynamics, src returns full state (including vels etc.)
%   -added "flag," which marks when a restart occurs
% Revised: 02/16/12
%   -heavy revision: functionalized, rationalized, etc.
% Revised: 02/03/11
%   -added varargin for biases and gains on all modalities
%   -made the appropriate other changes to accomodate these
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme
% Revised: 12/06/10 (happy b'day)
%   -changed to accomodate variable-dimensional stimuli
% Revised: 11/29/10
%   -changed to accomodate x in *true* coordinates
% Revised: 08/23/10
%   -added varargin for funny params
% Revised: 07/05/10
%   -consolidated params into setParams
%   -eliminated randpos fxn :-|
% Revised: 06/07/10
%   -accounted for edge effects, normalized the grid, etc.
% Created: 06/04/10
%   by JGM
%-------------------------------------------------------------------------%

% get "non-standard" variables
[Q,params] = specialvars(varargin,params);

% check for decoupling
Q.DECOUPLED = rand(Nbatches*params.Ncases,1) >= params.crln;

% say
if Q.VERBOSE; describeData(Nbatches,Q,params); end

% generate source vars
[S,Q] = getStimuli(Q,Nbatches,params);
    
% compute the non-source var and encode all in GTPNs
[R,G] = encodeStimuli(S,Q,Nbatches,params);

% matrix -> tensor
[R,G] = shortdata(params.Ncases,3,R,G);
S = shortdata(params.Ncases,4,S);
varargout{1} = R;
varargout{2} = S;
varargout{3} = G;
varargout{4} = Q;



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function describeData(Nbatches,Q,params)

fprintf('\n\n');
fprintf('Generating %iD data in %i modalities:\n',params.Ndims,params.Nmods);
fprintf('   modality correlation: %1.2f\n',params.crln);
fprintf('   variable gain:        +/- %.1f%%\n',100*params.swing);
fprintf('   tuning-curve max.:    %i spikes/unit time\n',Q.xpctG);
fprintf('   neurons/dimension:    %i\n',params.N);
fprintf('   number of data:       %i cases/batch for %i batches = %i ',...
    params.Ncases,Nbatches,params.Ncases*Nbatches);
fprintf('training vectors\n');
fprintf('\n\n');
%03.f
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Q,params] = specialvars(audibles,params)

% params
Ndims = params.Ndims;
Nmods = params.Nmods;
if strcmp(params.machine,'domestica')
    yrclass = 'gpuArray';
else
    yrclass = 'double';
end

% First search through params for special variables...
% (1) source modality (default = prop)
Q.SOURCE = 'Joint-Angle';
if isfield(params,'NS'), Q.SOURCE = params.NS; end
% (2) source sampling (default = uniform)
Q.srcSmplFxn = @(M)(scalefxn(rand(Ndims,M,yrclass),...
    zeros(Ndims,1,yrclass),ones(Ndims,1,yrclass),...
    params.smin(:,strcmp(params.mods,params.NS)),...
    params.smax(:,strcmp(params.mods,params.NS)))');
% if strcmp(Q.SOURCE,'Hand-Position')&&(Ndims==2)
%     Q.srcSmplFxn = @(M)(vissampler(M,params));
% end

if isfield(params,'p0')
    p0mu = params.p0.mu; p0cov = params.p0.cov;
    if det(p0cov) == 0; p0STD = 0; else p0STD = sqrtm(p0cov); end
    Q.srcSmplFxn = @(M)(bsxfun(@plus,...
        gpuArray(p0STD)*randn(params.Ndims,M),p0mu,yrclass)');
    %%% (p0STD*randn(params.Ndims,M) + repmat(p0mu,1,M))');
end

% (3) gain means ("defauts" = params.g)
Q.xpctG = params.g*ones(1,Nmods,yrclass);
if isfield(params,'gains'); Q.xpctG = params.gains; end


% set defaults
params.crln = 1;
Q.VERBOSE = 1;
Q.biases = zeros(Ndims,Nmods,yrclass);


% ...now seach through the "audibles" for special varibles
for i = 1:2:length(audibles)
    switch audibles{i}
        case 'swing'
            params.swing = audibles{i+1};
        case 'correlation'
            params.crln = audibles{i+1};
        case 'visbias'
            Q.biases(:,strcmp(params.mods,'Hand-Position')) = audibles{i+1};
        case 'propbias'
            Q.biases(:,strcmp(params.mods,'Joint-Angle')) = audibles{i+1};
        case 'eyebias'
            Q.biases(:,strcmp(params.mods,'Gaze-Angle')) = audibles{i+1};
        case 'visgain'
            Q.xpctG(strcmp(params.mods,'Hand-Position')) = audibles{i+1};
        case 'propgain'
            Q.xpctG(strcmp(params.mods,'Joint-Angle')) = audibles{i+1};
        case 'eyegain'
            Q.xpctG(strcmp(params.mods,'Gaze-Angle')) = audibles{i+1};
        case 'prior'
            p0mu = audibles{i+1}.mu(:);     % overwrites the prior in
            p0cov = audibles{i+1}.cov;      %  params, if there is one
            if det(p0cov) == 0; p0STD = 0; else p0STD = sqrtm(p0cov); end
            Q.srcSmplFxn = @(M)(bsxfun(@plus,p0STD*randn(params.Ndims,M,yrclass),p0mu)');
            %%% @(M)((p0STD*randn(Ndims,M,yrclass) +repmat(p0mu,1,M))');
        case 'source'
            Q.SOURCE = audibles{i+1};
        case 'verbose'
            Q.VERBOSE = audibles{i+1};
        case 'prevparams'
            Q.paramsOLD = audibles{i+1};
        case 'prevweights'
            Q.wtsOLD = audibles{i+1};
        case 'dbndepth'
            Q.iRBM = audibles{i+1};
        case 'dbnwts'
            Q.wts = audibles{i+1};
        case 'stimuli'
            Q.stims = audibles{i+1};
        case 'deadneurons'
            Q.deadInds = audibles{i+1};
        case 'gains'
            Q.gains = audibles{i+1};
        case 'ICs'
            Q.s0 = audibles{i+1};
        case 'RNNwts'
            Q.RNNwts = audibles{i+1};
        otherwise
            fprintf('\n unrecognized option for DATAGENPP!!\n');
    end
end

Q.xpctG = Q.xpctG*ones(1,1,yrclass);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [S,Q] = getStimuli(Q,Nbatches,params)

% init
Ndims = params.Ndims;
Nexamples = Nbatches*params.Ncases;
uniformSmplFxn = @(M,smin,smax)(...
    scalefxn(rand(Ndims,M,'like',Q.xpctG),zeros(Ndims,1,'like',Q.xpctG),...
    ones(Ndims,1,'like',Q.xpctG),smin,smax)');

if isfield(Q,'stims')
    S = Q.stims; 
else

    % different models call for stimuli
    switch params.MODEL
        
        case {'2Dinteg','2DintegDecoupled'}
            switch Q.SOURCE
                case 'Hand-Position'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.thmin,params.thmax));
                    convertfxn = @(S)(IK2link(S,params,1));
                    [vis,prop] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                case 'Joint-Angle'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.posmin,params.posmax));
                    convertfxn = @(S)(FK2link(S,params,1));
                    [prop,vis] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                otherwise
                    error('strange modality for this model -- jgm');
            end
            S(:,:,strcmp(params.mods,'Hand-Position')) = vis;
            S(:,:,strcmp(params.mods,'Joint-Angle')) = prop;
            
        
        case {'2Dtwoarms'}
            nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.smin(:,1),params.smax(:,1)));
            convertfxn = @(S)(S);
            [propL,propR] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
            
            S(:,:,strcmp(params.mods,'Joint-Angle-Left')) = propL;
            S(:,:,strcmp(params.mods,'Joint-Angle-Right')) = propR;
        
            
        case '1Dinteg'
            switch Q.SOURCE
                case 'Hand-Position'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.thmin,params.thmax));
                    convertfxn = @(S)(IK2link(S,params,1));
                    [vis,prop] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                case 'Joint-Angle'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.posmin,params.posmax));
                    convertfxn = @(S)(FK2link(S,params,1));
                    [prop,vis] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                otherwise
                    error('strange modality for this model -- jgm');
            end
            S(:,:,strcmp(params.mods,'Hand-Position')) = vis;
            S(:,:,strcmp(params.mods,'Joint-Angle')) = prop;
            
            
            
        case '1Daddition'
            eye = uniformSmplFxn(Nexamples,params.eyemin,params.eyemax);
            switch Q.SOURCE
                case 'Hand-Position'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.thmin,params.thmax));
                    convertfxn = @(S)(IK2link(S-eye,params,1));
                    [vis,prop] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                case 'Joint-Angle'
                    nonsrcSmplFxn = uniformSmplFxn(N,params.posmin,params.posmax);
                    convertfxn = @(S)(FK2link(S,params,1)+eye);
                    [prop,vis] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                otherwise
                    error('strange modality for this model -- jgm');
            end
            S(:,:,strcmp(params.mods,'Hand-Position')) = vis;
            S(:,:,strcmp(params.mods,'Joint-Angle')) = prop;
            S(:,:,strcmp(params.mods,'Gaze-Angle')) = eye;
            
            
            
        case '2Daddition'
            eye = uniformSmplFxn(Nexamples,params.eyemin,params.eyemax);
            switch Q.SOURCE
                case 'Hand-Position'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.thmin,params.thmax));
                    convertfxn = @(S)(IK2link(S-eye,params,1));
                    [vis,prop] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                case 'Joint-Angle'
                    nonsrcSmplFxn = @(N)(uniformSmplFxn(N,params.posmin,params.posmax));
                    convertfxn = @(S)(FK2link(S,params,1)+eye);
                    [prop,vis] = getStimuliCore(Q.srcSmplFxn,nonsrcSmplFxn,convertfxn,Q.DECOUPLED);
                otherwise
                    error('strange modality for this model -- jgm');
            end
            S(:,:,strcmp(params.mods,'Hand-Position')) = vis;
            S(:,:,strcmp(params.mods,'Joint-Angle')) = prop;
            S(:,:,strcmp(params.mods,'Gaze-Angle')) = eye;
          
            
            
        case 'MCDexperiment'
            
            dmin = params.experiment.dmin;
            dmax = params.experiment.dmax;
            r = sqrt(scalefxn(rand(Nexamples,1,'like',Q.xpctG),0,1,dmin^2,dmax^2));
            th = scalefxn(rand(Nexamples,1,'like',Q.xpctG),0,1,-pi,pi);

            
            %%% replace with pol2cart ??
            S(:,1,strcmp(params.mods,'Motion-Dots')) = r.*cos(th);
            S(:,2,strcmp(params.mods,'Motion-Dots')) = r.*sin(th);
            S(:,:,strcmp(params.mods,'ICMS')) = [th, r];
            
            
        case 'MCDdarpa'
            
            dmin = params.experiment.dmin;
            dmax = params.experiment.dmax;
            r = sqrt(scalefxn(rand(Nexamples,1,'like',Q.xpctG),0,1,dmin^2,dmax^2));
            th = scalefxn(rand(Nexamples,1,'like',Q.xpctG),0,1,-pi,pi);

            
            %%% replace with pol2cart ??
            S(:,1,strcmp(params.mods,'Motion-Dots')) = r.*cos(th);
            S(:,2,strcmp(params.mods,'Motion-Dots')) = r.*sin(th);
            S(:,1,strcmp(params.mods,'ICMS')) = r.*cos(th);
            S(:,2,strcmp(params.mods,'ICMS')) = r.*sin(th);
            
                       
        case 'HierL2'
            [D0,S0,G0] = DATAGENPP(Nbatches,Q.paramsOLD);
            Q.D = D0;
            Q.S = S0;
            Q.G = G0;
            
            %%% this could stand to be more general....
            switch params.typeUnits{1}
                case 'PB'
                    newinputind = 1;
                case 'BP'
                    newinputind = 2;
                otherwise
                    error('unexpected unit types -- jgm');
            end
            S = longdata(S0);
            S = S(:,:,newinputind); 
            
        case {'1DrEFH','2DrEFH','1DtRBM','1DrEFHwithEC','2DrEFHwithEC',...
                '3DrEFH','1DrEFHbern','1DrEFHwithECdelayed',...
                'HHSreachData','rEFH','TRBM','RTRBM'}
            
            if isfield(Q,'s0')
                params.dynamics.muX0 = Q.s0(:,1:end/2);
                params.dynamics.SigmaX0 = zeros(size(params.dynamics.SigmaX0));
                params.dynamics.muV0 = Q.s0(:,(end/2+1:end));
                params.dynamics.SigmaV0 = zeros(size(params.dynamics.SigmaV0));
            end
            [S,filterdata] = trajectoryGen(Nbatches,params);
            Q.filterdata = filterdata;
            
            
        case 'BMMdata'
            %%%% you don't really need to load anything here, although S has to
            %%%% be assigned something....
            
        otherwise
            error('unrecognized model -- jgm');
            
    end
    
    
    % check for out-of-range values of the stimulus
    %%%%% put this back in!  It's just that, on the torus, you think you
    %%%%% can get values slightly over the max......
%     for iMod = 1:size(S,3)
%         for iExample = 1:size(S,1)
%             s = squeeze(S(iExample,:,iMod));
%             outofrange(s,params.smin(:,iMod),params.smax(:,iMod),...
%                 params.mods{iMod},1,params);
%         end
%     end

end




end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [ssrc,snonsrc] = getStimuliCore(srcSmplFxn,nonsrcSmplFxn,...
    convertfxn,DECOUPLED)

% source space
ssrc = srcSmplFxn(length(DECOUPLED)); 

% non-source space
snonsrc(~DECOUPLED,:) = convertfxn(ssrc(~DECOUPLED,:));
if any(DECOUPLED), snonsrc(DECOUPLED,:) = nonsrcSmplFxn(sum(DECOUPLED)); end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [R,G] = encodeStimuli(S,Q,Nbatches,params)
% The master function for all the different ways we might encode S into
% neural responses D.

%%%%%%
% Still an open question whether you should count the feedback as a
% "modality" (Nmods)--right now you don't.
%%%%%%


% different models call for different data
switch params.MODEL
    
    
    case {'2Dinteg','1Daddition','1Dinteg','2Daddition',...
            '2DintegDecoupled','MCDdarpa','2Dtwoarms'}
    
        [R,G] = stdStimulusEncoding(S,Q,Nbatches,params.typeUnits{1},params);
        
        
    case 'HierL2'
        
        % the third input
        Nvis = params.numsUnits(1);
        switch params.typeUnits{1}
            case 'PB', 
                FXN = 'Poisson'; 
                inds3 = 1:(Nvis-params.t);
                indsHids = (params.t+1):Nvis;
            case 'BP', 
                FXN = 'Poisson'; 
                inds3 = (params.t+1):Nvis;
                indsHids = 1:(Nvis-params.t);
            otherwise, error('unexpected unit types in hierarchy -- jgm');
        end
        [R3,G] = stdStimulusEncoding(S,Q,Nbatches,FXN,params);
        
        
        % the first and second inputs
        R1and2 = longdata(Q.D);
        wts0 = Q.wtsOLD;
        params0 = Q.paramsOLD;
        smpls = params0.smpls;
        
        % generate the half of second level data given by the first level
        means = feedforward(R1and2,wts0{1}(1:end-1,:),wts0{1}(end,:),...
            params0.typeUnits{2},params0);
        states = sampler(means,params0.typeUnits{2},params0);
        for i = 1:smpls-1
            states = states + sampler(means,params0.typeUnits{2},params0);
        end
        Rhiddens = states/smpls;
        
        % concatenate
        R(:,inds3) = R3; R(:,indsHids) = Rhiddens;
        
        
        
    case 'BMMdata'
        %%%% R,G = load some data....
        
        
        
    case {'1DrEFH','2DrEFH','1DtRBM','3DrEFH','1DrEFHwithEC',...
            '2DrEFHwithEC','1DrEFHbern','1DrEFHwithECdelayed',...
            'HHSreachData','rEFH','TRBM','RTRBM'}
        
        % what kind of input units?
        switch params.typeUnits{1}
            case 'PB', FXN = 'Poisson';
            case 'BP', FXN = 'Poisson';
            case 'BG', FXN = 'Gaussian';
            case 'Bernoulli', FXN = 'Bernoulli';
            otherwise, error('unexpected unit types in rEFH -- jgm');
        end
        
        % the sensory data
        if strcmp(params.dynamics.walls,'wrapping')
            [Rsensory,G] =...
                encodeToroidalStimuli(S,Q,Nbatches,FXN,params);
        else
            [Rsensory,G] =...
                stdStimulusEncoding(S,Q,Nbatches,FXN,params);
        end
        
        Nhids = params.numsUnits(2);
        %%% Nexamples = Nbatches*params.Ncases;
        Ncases = params.Ncases;
        %%% Nbatches = Nexamples/Ncases;  what was this for?
        
        
        
        % create the recurrent data---perhaps "blank"
        if isfield(Q,'RNNwts')
            RrecurInit = zeros(Ncases,Nhids,'like',Rsensory);
            Rrecur = RNNforwardpass(...
                shortdata(Ncases,3,Rsensory),RrecurInit,...
                Q.RNNwts(1:end-1,:),Q.RNNwts(end,:),...
                params.typeUnits{2},params);
            %%% params.typeUnits{2} should be fixed for DBNs...
        else
            %%% a little hacky....
            Rrecur = -ones(Ncases,Nhids,Nbatches); %%% 'like' ???
            if ~isfield(Q,'s0')
                Rrecur(:,:,1) = 0;
                fprintf('setting initial "recurrence" to a vector of zeros\n\n');
            end
        end
        
        % write special information into time steps greater than 1
        if strcmp(params.dynamics.walls,'resetting')
            for iBatch = 2:Nbatches
                Rrecur(filterdata(iBatch).RESTART,:,iBatch) = 0;
            end
        end
        
        % concatenate
        Rrecur = longdata(Rrecur);
        switch params.typeUnits{1}
            case 'PB'
                R = cat(2,Rsensory,Rrecur);
            case {'BP','GP','BG','Bernoulli'}
                %%% so Bern-Bern rEFH also puts the recurrent on the left!
                R = cat(2,Rrecur,Rsensory);
            otherwise
                error('unexpected unit types in recurrent EFH -- jgm');
        end
        
        
        
    case 'MCDexperiment'
        
        Sp = S(:,:,strcmp(params.mods,'ICMS'));
        Sc = S(:,:,strcmp(params.mods,'Motion-Dots'));
        F = params.experiment.electrodeFunc(Sp(:,1),Sp(:,2));
        Nns = params.experiment.neuronsPerElectrode;
        Ricms = sampler(repmat(F,[1,1,Nns]),params.typeUnits{1},params);
        Ricms = reshape(Ricms,[size(F,1),size(F,2)*Nns]);
        Gicms = params.g*ones(size(Ricms,1),1,'like',Ricms);

        params.Nmods = 1;
        [Rvis, Gvis] = stdStimulusEncoding(Sc,Q,Nbatches,...
            params.typeUnits{1},params);
        
        R = cat(2,Rvis,Ricms);
        G = cat(2,Gvis,Gicms);
        %%% NB: this assumes that VIS is on the left!!
        
        
    otherwise
        error('unrecognized model -- jgm');
end



if isfield(Q,'iRBM')    
    means = R;
    wts = Q.wts;
    
    for layer = 1:(Q.iRBM-1)
        HIDFXN = params.typeUnits{layer+1};
        means = feedforward(means,wts{layer}(1:end-1,:),...
            wts{layer}(end,:),HIDFXN,params);
%         if (layer==1)&&isfield(params,'dynamics')&&(iBatch<Nbatches)
%             %%% data for the first layer of RBM number 2, for dRBM
%             Din(:,1:numsUnits(2),iBatch+1) = means;
%         end
    end
end

% kill certain neurons
if isfield(Q,'deadInds'), R(:,Q.deadInds) = 0; end



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [R,G] = stdStimulusEncoding(S,Q,Nbatches,FXN,params)


% init
Nmods = params.Nmods;
Ncases = params.Ncases;
Ndims = params.Ndims;
N = params.N;
smin = params.smin;
smax = params.smax;
patchmin = params.margin*ones(1,Ndims);
patchmax = patchmin + params.respLength;
Nexamples = Nbatches*Ncases;

% draw gains
if isfield(Q,'gains')
    G = Q.gains; 
else
    G = unifSmplAboutCenter(Q.xpctG,params.swing,Nexamples);
end

% draw population responses
R = zeros(Nexamples,Nmods*N^Ndims,'like',S);
for iMod = 1:Nmods
    
    % prepare
    b = Q.biases(:,iMod);
    thissmin = smin(:,iMod);    thissmax = smax(:,iMod);
    
    theseS = bsxfun(@plus,S(:,:,iMod)',b);
    sPatch = scalefxn(theseS,thissmin,thissmax,patchmin,patchmax)';
    inds = ((iMod-1)*N^Ndims + 1):(iMod*N^Ndims);
    
    % encode the stimuli in GTPNs
    R(:,inds) = PPCencode(sPatch,G(:,iMod),FXN,params);
end


% mark the decoupled vectors
R(:,params.N) = Q.DECOUPLED*params.g + ~Q.DECOUPLED.*R(:,params.N);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [R,G] = encodeToroidalStimuli(S,Q,Nbatches,FXN,params)
% recall that S is 40,000 x Ndims x Nmods; smin, smax are Ndims x Nmods

% Ns
Nexamples = size(S,1);
Ndims = params.Ndims;
Nmods = params.Nmods;
N = params.N;

% first generate all the gains, which must be fixed across each recursion
if isfield(Q,'gains')
    G = Q.gains; 
else
    G = unifSmplAboutCenter(Q.xpctG,params.swing,Nexamples);
    Q.gains = G;
end
%%% NB: Q.gains is what is used in stdStimulusEncoding; so it is gauranteed
%%% to be the same all across all nine (or whatever) repeats of the
%%% toroidal encoding.

% now wrap the stimuli into encoding space
smin = reshape(params.smin,[1,Ndims,Nmods]);    % reshape these for bsxfun
smax = reshape(params.smax,[1,Ndims,Nmods]);    % 
srange = N/(N-1)*(smax - smin);                 % confusing, to be sure
S = bsxfun(@plus,bsxfun(@mod,bsxfun(@minus,S,smin),srange),smin);

% then displace them by length of range(s) for the recursive thing
S = bsxfun(@minus,S,srange);

% since all the gains have been sampled, don't resample
[R,~] = encodeRecursively(S,Q,Nbatches,params,1,FXN);

end
%-------------------------------------------------------------------------%
 
 
%-------------------------------------------------------------------------%
function [R,S] = encodeRecursively(S,Q,Nbatches,params,dim,FXN)

% init
R = 0;
Ndims = params.Ndims;
Nmods = params.Nmods;
N = params.N;

% (re-)calculate srangemat
srange = reshape(N/(N-1)*(params.smax - params.smin),[1,Ndims,Nmods]);

if dim <= Ndims
    thisS = S;
    for iWrap = 1:3
        [thisR,thisS] = encodeRecursively(thisS,Q,Nbatches,params,dim+1,FXN);
        R = R + thisR;
        
        thisS(:,dim,:) = bsxfun(@plus,thisS(:,dim,:),srange(1,dim,:));
    end
else
    [R,~] = stdStimulusEncoding(S,Q,Nbatches,FXN,params);
end
    

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function inputUnitType = setInputUnitType(typeUnits,Nmods)

switch typeUnits
    case 'BP'
        inputUnitType{1} = 'Bernoulli';
        inputUnitType{2} = 'Poisson';
    case 'PB'
        inputUnitType{1} = 'Poisson';
        inputUnitType{2} = 'Bernoulli';
    otherwise
        for i = 1:Nmods
            inputUnitType{i} = typeUnits;
        end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function smpl = priorsmplr(Ndims,center,numstd)

if numstd
    onestddev = eye(Ndims)/numstd;              % whole range w/i numstd std
    smpl = onestddev*randn(Ndims,1) + center;
else
    smpl = center;
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Ncases = params.Ncases;
% Tstretch = 10;
% gmax = (1 + params.swing)*params.g;
% gmin = (1 - params.swing)*params.g;     
% gainSchedule = ones(Tstretch,1)*repmat([gmin gmax],1,Nbatches/Tstretch/2);
% g = reshape(repmat(gainSchedule(:)',[Ncases,1]),[Ncases,1,Nbatches]);
% G = longdata(g);
%
% [R,S] = DATAGENPP(1000,params,'gains',G);
%-------------------------------------------------------------------------%
