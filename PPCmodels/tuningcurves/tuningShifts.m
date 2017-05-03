function tuningShifts(shftstd,wts,params)
% tuningShifts(shftstd,wts,params)
%   tuningShifts plots the tuning of (some of) the hidden units for both
%   shifted an unshifted data.  Do input shifts merely shift the hidden
%   layer tuning curves, or do they change their shapes as well?  Look and
%   see!
%
% MRF requested to see these results for his paper
%
% first load:
%   load results/numhidswts/Std050.mat
% then run:
%   tuningShifts(shftstd,wts,params);

%-------------------------------------------------------------------------%
% Created: 08/02/12
%   by JGM
%-------------------------------------------------------------------------%

% shftstd = 5;

%%%%%
% Currently (1/1/17) broken, but should be straightforward to fix...
%%%%% 


% init
M = 40;
Ndims = params.Ndims;
Nmods = length(params.mods);
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};
gains = mean([params.gmin; params.gmax]);
params.gmin = gains;
params.gmax = gains;
indices = reshape(1:Nmods*Ndims,Ndims,Nmods);
indout = indices(:,strcmp(params.mods,params.NS));
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

% get "baseline" tuning curves
[X,Q] = getStimuliTiled(M^2,dataclass,params);
R = params.getData(X,Q);
V = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,hidNums,params);

% transform shftstd (of "input covariance") into units of params.NS
shft = std2shft(shftstd,X,R,params);

% get "shift-data" tuning curves
Xshft = X;
Xshft(:,indout) = X(:,indout) + repmat(shft',size(X,1),1);
Q.G = repmat(mean([params.gmin; params.gmax]),[M^2,1]);
R = params.getData(Xshfit,Q);
Vshft = invParamMap(Rshft,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,...
    hidNums,params);

% get the integrated estimate from the input populations
xopt = getOptimalStims(Xshft,R,params);


% plot tuning curves
plotTCs(X(:,3:4),V,params);
% plotTCs(X(:,3:4),cat(3,V,Vshft),params);
% plotTCs(Xshft(:,3:4),Vshft,params);
plotTCs(xopt,Vshft,params);





end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function shft = std2shft(shftstd,X,R,params)

keyboard
% compute input error covariance, CovV + CovP (in neutral space) for shft
SINSMerr = covInCalc(R,X,params);
ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;

%%% hard-coded direction!
nDirections = 8;
iDirection = 3;

% compute shift (in neutral space)
phi = 2*pi*iDirection/nDirections;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
shft = sqrtm(ErrCovIn)*R*[shftstd; 0];



end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function xopt = getOptimalStims(Xshft,R,params)
%%% this assumes length(params.mods) = 2!


% init
Nexamples = size(R,1);
Ndims = params.Ndims;

% loop
xopt = zeros(Nexamples,Ndims);
for iExample = 1:Nexamples
    
    rr = R(iExample,:);
    xx = Xshft(iExample,:);
    
    % get estimates and covariances
    shatL = decoder(rr,params);
    shatN = estGather(shatL,params);
    SINSMerr = covInCalc(rr,xx,params);
    
    % combine optimally
    pcnV = SINSMerr{1}.cov;
    pcnP = SINSMerr{2}.cov;
    WV = (pcnV+pcnP)\pcnV;
    WP = (pcnV+pcnP)\pcnP;
    
    % x = estGather(reshape(Xshft(iExample,:),Ndims,Nmods),params);
    xopt(iExample,:) = [WV WP]*shatN(:);
end

end
%-------------------------------------------------------------------------%