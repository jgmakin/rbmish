function varargout = optMeanCheckPP(D,params,MODE,varargin)
%%% vopt = optMeanCheckPP(D,xtrue,wts,params,s,MODE)
%
% OPTMEANCHECKPP    Compare optimal and RBMed means, trial by trial
%   OPTMEANCHECKPP computes the trial-by-trial estimate of the visual
%   stimulus given by (1) the rbm (vRBM) and (2) the optimal procedure
%   (vopt).  The latter can be computed using any of three covarince
%   matrices: the actual error covariances (from estimatorStats and
%   contained in s), the ideal error covariances (computed via the Fisher
%   info), and the trial-by-trial covariance (computed via the *observed*
%   per trial Fisher info): MODE = 'errorCov'/'expectedFI'/'observedFI'.
%
%   EXAMPLES:
%       [vopt vRBM x] = optMeanCheckPP(D0,params,'expectedFI',x0,wts);
%       [vopt vRBM x] = optMeanCheckPP(D0,params,'observedFI',x0,wts);
%       [vopt vRBM x] = optMeanCheckPP(D0,params,'errorCov',x0,wts,sV0);
%       vopt = optMeanCheckPP(D0,params,'expectedFI');
%       vopt = optMeanCheckPP(D0,params,'observedFI');
%       vopt = optMeanCheckPP(D0,params,'errorCov',sV0);
%
%   NB this appears to be broken

%-------------------------------------------------------------------------%
% Revised: 12/12/13
%   -changed to accommodate new orientation of Dlong and Xlong....
% Created: 11/30/10
%   by JGM
%-------------------------------------------------------------------------%

% extract params
L1 = params.L1;  L2 = params.L2;
pmin = params.thmin;
pmax = params.thmax;
vmin = params.posmin;
vmax = params.posmax;

% parloop init
D0 = longdata(D);
Nexamples = size(D0,1);

if nargin > 3
    if nargin == 4
        if ~strcmp(MODE,'errorCov')
            fprintf('one extra input argument but not in errorCov mode');
            error('(expected sV) --- jgm');
        else
            s = varargin{3};
        end
    else
        xtrue = varargin{1};
        wts = varargin{2};
        if nargin == 6
            s = varargin{3};
        end
        x = reshape(shiftdim(xtrue,1),size(xtrue,2),Nexamples)';
        tic
        [avgErr,D1] = updown(D0,wts,params,'means');
        toc
        vRBM = zeros(2,Nexamples);
        varargout{3} = x;
    end
end
%     if nargout ~= 3
%         error('wrong number of outputs for these inputs -- jgm');
%     else

switch MODE
    case 'errorCov'
        FIvv = inv(s.covVi);
        FIpp = inv(s.covPi);
    case 'expectedFI'
        FIvv = PPCexpectedFI([vmin,vmax],params.g,params);
        FIpp = PPCexpectedFI([pmin,pmax],params.g,params);
    case 'observedFI'
        FIvv = 0;
        FIpp = 0;
    otherwise
        error('unrecognized mode -- jgm');
end

% parloop
vopt = zeros(2,Nexamples);
HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED, pool = parpool(4); end
parfor iExample = 1:Nexamples
    
    % compute pre-updated estimates
    T0 = displayshape(D0(iExample,:),params);
    v0 = decode(T0{1},[vmin vmax],params,{'CoM'});
    p0 = decode(T0{2},[pmin pmax],params,{'CoM'});
    
    % compute optimal x estimate
    J =[-L1*cos(p0(1)) - L2*cos(p0(1)+p0(2)), -L2*cos(p0(1)+p0(2));...
        -L1*sin(p0(1)) - L2*sin(p0(1)+p0(2)), -L2*sin(p0(1)+p0(2))];
    switch MODE
        case 'observedFI'
            FIovv = PPCobservedFI(V0,v0,[vmin,vmax],params);
            FIopv = J*PPCobservedFI(P0,p0,[pmin,pmax],params)*J';
            vopt(:,iExample) = (FIovv + FIopv)\FIovv*v0 +...
                        (FIovv + FIopv)\FIopv*FK2link(p0,params,0)';
        case {'errorCov','expectedFI'}
            FIpv = J*FIpp*J';
            vopt(:,iExample) = (FIvv + FIpv)\FIvv*v0 +...
                        (FIvv + FIpv)\FIpv*FK2link(p0,params,0)';
        otherwise
            error('unrecognized mode -- jgm');
    end
    
    if nargin > 0
        % compute model's x estimate
        T1 = displayshape(D1(iExample,:),params);
        vRBM(:,iExample) = decode(T1{1},[vmin vmax],params,{'CoM'});
    end
    
end
if HADBEENCLOSED, delete(pool); end


varargout{1} = vopt;
if ~isempty(varargin)
    varargout{2} = vRBM;
end

end







