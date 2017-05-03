function [S,Q] = getStimuliTiled(Nexamples,yrclass,params,varargin)
%%%function S = getUniformStims(M,datatype,mods,roboparams,varargin)
% Creates a uniform MxM grid of joint angles, and the corresponding
% end-effector positions, in the usual format.

%-------------------------------------------------------------------------%
% Revised: 01/01/17
%   -renamed from getUniformStims -> getStimuliTiled
%   -heavily modified to parallel getStimuliMultisensoryIntegration.  To
%   process is not complete, and this fxn is left rather inelegant.
% Revised: 12/16/13
%   -X -> S and corresponding changes
% Created: 07/27/12
%   by JGM
%-------------------------------------------------------------------------%


% params
mods = params.mods;
roboparams = params.roboparams;
Ndims = length(roboparams.thmin);

% some useful vectors
M = sqrt(Nexamples);
if round(M)~=M
    error('getStimuliUniform requires Nexamples be a perfect square\n');
end
z = linspace(zeros(1,yrclass),ones(1,yrclass),M);
oo = ones(Ndims,1,yrclass);

if any(strcmp(mods,'Gaze-Angle'))
    
    % init
    posmin = roboparams.posmin;     posmax = roboparams.posmax;
    eyemin = roboparams.eyemin;     eyemax = roboparams.eyemax;
    eyespacefraction = defaulter('eyefraction',1/3,varargin{:});
    emin = eyemin*eyespacefraction;
    emax = eyemax*eyespacefraction;
    
    % Tbody increases left to right
    Tbody = permute(repmat(linspace(posmin+emax,posmax+emin,M)',...
        [1 1 1 4]),[4,2,3,1]);
    
    % *gaze* increases from *bottom to top* (eye *decreases*)
    gaze = repmat(-linspace(emin,emax,M)',[1 1 1 4]);
    
    % VIS, PROP, EYE
    S = cat(3,Tbody - gaze,IK2link(Tbody,roboparams,1),-gaze);
    S = longdata(S);
else
    
    % init
    thmin = roboparams.thmin;
    thmax = roboparams.thmax;
    
    alllattices = ndlattice(Ndims,z);
    prop = scalefxn(alllattices,0*oo,1*oo,thmin,thmax);
    S(:,:,strcmp(mods,'Joint-Angle')) = prop;
    
    if isfield(mods,'Hand-Position')
        vis = FK2link(prop,roboparams,1);
        S(:,:,strcmp(mods,'Hand-Position')) = vis;
    end
    
end



% now populate Q
% gains
uniformSmplFxn = @(M,xmin,xmax)(UniformNormalDiracSampler(...
    zeros(length(xmin),1,yrclass),Inf,M,xmin,xmax,0));
Q.G = defaulter('gains',uniformSmplFxn(Nexamples,params.gmin,params.gmax),...
    varargin{:});

% biases (default to none)
Q.biases = zeros(size(params.smin),yrclass);
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'visbias'
            Q.biases(:,strcmp(mods,'Hand-Position')) = varargin{i+1};
        case 'propbias'
            Q.biases(:,strcmp(mods,'Joint-Angle')) = varargin{i+1};
        case 'eyebias'
            Q.biases(:,strcmp(mods,'Gaze-Angle')) = varargin{i+1};
    end
end

% "coupling"
Q.DECOUPLED = zeros(Nexamples,1,yrclass);

end

