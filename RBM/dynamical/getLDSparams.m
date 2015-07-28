function LDSparams = getLDSparams(params,learningtype,varargin)
% A few different methods for getting the params of a "linear dynamical
% system" (linear, time-invariant system).
%
% USAGE: LDSparams = getLDSparams(T,params,learningtype)
%
% NB: we're assuming throughout, wlog, that B = I and the transition noise
% has covariance SigmaX.
%
% NB: this file should take care of all the things that are peculiar to
% population coding.  The subfunctions should work for *any* LDS.
%
% NB: that this will fail if params doesn't have a field "dynamics"

%-------------------------------------------------------------------------%\
% Revised: 08/26/14
%   -changed case "true" to treat "efference copied" controls as states
% Revised: 01/10/14
%   -getXandY -> getLDSdata, and associated changes, esp. to accommodate
%   controlled systems
% Revised: 10/30/13
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%



switch learningtype
    case 'observed'

        if ~isempty(varargin)                               % get data
            LDSdata = varargin{1}; 
        else
            LDSdata = getLDSdata(params);                  
        end
        LDSparams = learnfullyobservedLDS(LDSdata);         % fit

        
    case 'true'
       
        LDSparams = getLDSparamsTrue(params);
        
    case 'nolearning'
        
        %%% currently messed up
        LDSparams = manualAdjustLDSparams(params);
        LDSparams.T = params.dynamics.T;
        %%%
        
    case 'EM'
        
        % this is the version that enforces noisly observed ctrls (zeros
        % out xpctTs)
        params.dynamics.Nstateobsvs = size(params.dynamics.C,1);
        Nstates = size(params.dynamics.A,2);        
        LDSparams = EM4LDS(Nstates,params);             % fit
        
    case 'presaved'
        
        switch params.MODEL
            case '1DrEFH'
                % load('dynamical/LDSparamsEM1Dwrapping140519.mat','LDSparamsEM');
                % optimized x100
                load('dynamical/LDSparamsEM1Dwrapping140520.mat','LDSparamsEM');
                %%% this is for no damping/higher velocity noise
            case '2DrEFH'
                load('dynamical/LDSparamsEM2Dwrapping140507.mat','LDSparamsEM');
                % optimized x100
                %%% No!! 140507 had a bug in the code
            case '1DrEFHwithEC'
                load('dynamical\LDSparamsEM1Dcontrolledfull140828.mat','LDSparamsEM');
                % optimized x100
            case '2DrEFHwithEC'
                load('dynamical/LDSparamsEM2Drepelling140422.mat','LDSparamsEM');
            otherwise
                error('strange number of dimensions! -- jgm');
        end
        LDSparams = LDSparamsEM;
        
    otherwise
        fprintf('whoops, not a recognized way to learn an LDS\n\n');
        LDSparams.error = 1;
end





end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function LDSparams = getLDSparamsTrue(params)
%%% TO DO:
% (1) if there's a visible control....
% (2) GPUify based on params.machine??
% (3) get rid of muX0, muV0 in params.dynamics?
%%%%%%%%%%%%%%



dynamics = params.dynamics;
%%%%% only works under the assumption of 2nd-order dynamics
Ndims = length(dynamics.muX0);
mu0 = [dynamics.muX0; dynamics.muV0];
Sigma0 = zeros(2*Ndims);
Sigma0(1:Ndims,1:Ndims) = dynamics.SigmaX0;
Sigma0(Ndims+1:end,Ndims+1:end) = dynamics.SigmaV0;


if any(strcmp(params.mods,'Efference-Copy'))
    
    % Ns
    [Ny,Nx] = size(dynamics.C);
    [Nv,Nu] = size(dynamics.H);
    Nstates = Nx+Nu;
    Nobsvs = Ny+Nv;
    
    % assemble the augmented state, augmented observations, etc.
    Gamma = zeros(Nstates,Nstates);
    Gamma(1:Nx,1:Nx) = dynamics.A;
    Gamma(1:Nx,(Nx+1):end) = dynamics.G;
    Lambda = zeros(Nobsvs,Nstates);
    Lambda(1:Ny,1:Nx) = dynamics.C;
    Lambda(Ny+1:end,Nx+1:end) = dynamics.H;
    nuZ = zeros(Nstates,1);
    try nuZ(1:Nx) = dynamics.muX; catch ME; end
    UpsilonZ = zeros(Nstates,Nstates);
    UpsilonZ(1:Nx,1:Nx) = dynamics.SigmaX;
    nu0 = zeros(Nstates,1);
    nu0(1:Nx) = mu0;
    Upsilon0 = zeros(Nstates,Nstates);
    Upsilon0(1:Nx,1:Nx) = Sigma0;
    
    % control cumulants, etc.
    if strcmp(params.dynamics.walls,'controlled')
        
        % then controls' cumulants are to be acquired empiricaly
        LDSdata = getLDSdata(params);
        U = longdata(LDSdata.Z(:,Nx+1:end,:));
        muU = mean(U);                  SigmaU = cov(U);
        muU0 = muU;                     SigmaU0 = SigmaU;
    else
        
        % the ctrl has real dynamics, whose params are in params.dynamics
        muU = zeros(Nu,1);              SigmaU = params.dynamics.SigmaU;
        muU0 = params.dynamics.muU0;    SigmaU0 = params.dynamics.SigmaU0;
        
        % fix state-transition matrix to account for dynamics in ctrl
        Gamma((Nx+1):end,(Nx+1):end) = params.dynamics.F;
    end
    nuZ(Nx+1:end) = muU;            UpsilonZ(Nx+1:end,Nx+1:end) = SigmaU;
    nu0(Nx+1:end) = muU0;           Upsilon0(Nx+1:end,Nx+1:end) = SigmaU0;
    
    
    % collect
    LDSparams.A = Gamma;
    LDSparams.C = Lambda;
    LDSparams.muX = nuZ;
    LDSparams.SigmaX = UpsilonZ;
    LDSparams.mu0 = nu0;
    LDSparams.Info0 = inv(Upsilon0);
    
else
    
    
    LDSparams.A = dynamics.A;
    LDSparams.C = dynamics.C;
    try LDSparams.muX = dynamics.muX;
    catch ME, LDSparams.muX = zeros(size(mu0));
    end
    LDSparams.SigmaX = dynamics.SigmaX;
    LDSparams.mu0 = mu0;
    LDSparams.Info0 = inv(Sigma0);
    
end
LDSparams.T = params.dynamics.T;
    
end
%-------------------------------------------------------------------------%
        
        
%-------------------------------------------------------------------------%
function LDSparams = manualAdjustLDSparams(params)
% this was your original method of "fixing" the KF parameters in the case
% of repelling walls.  The adjustment numbers were computed empirically.


% give away all the parameters
LDSparams.A = params.dynamics.A;
LDSparams.C = params.dynamics.C;
LDSparams.SigmaX = params.dynamics.SigmaX;
LDSparams.mu0 = [params.dynamics.muX0; params.dynamics.muV0];
LDSparams.Info0 = setInfoMatrix(params.dynamics.SigmaX0,...
    params.dynamics.SigmaV0,params.Ndims);


% then adjust the state-transition matrix...
A = LDSparams.A;
A(3,3) = 1-0.001;
A(4,4) = 1-0.001;
LDSparams.A = A;


% ...and the state-transition noise
SigmaX = LDSparams.SigmaX;
SigmaX(3:4,3:4) = SigmaX(3:4,3:4) + 1e-3*[0.1230 0.0098; 0.0098 0.2163];
%%% hacky: should invoke params.Ndims
LDSparams.SigmaX = SigmaX;



end
%-------------------------------------------------------------------------%
