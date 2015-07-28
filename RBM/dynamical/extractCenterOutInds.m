function cntrOutInds = extractCenterOutInds(S,params,TOPLOT)
% extractCenterOutInds
% 
% USAGE: 
%   cntrOutInds = extractCenterOutInds(LDSdataTest.S,params,0);
%
% This function was designed to work with params.MODEL = HHSreachData.  It
% takes the entire stimulus history, S, (Ncases x Ndims x Nmods x T), and
% extracts the starting and ending indices of only the center-out reaches;
% i.e., leaving out the out-center reaches.  Since these are not guaranteed
% to be of the same length, they are stored in cell arrays: cntrOutInds has
% length Ncases, each cell of which holds 2 x a variable number of indices
% ("2" b/c they're the start and end indices).

%-------------------------------------------------------------------------%
% Created: 03/18/15
%   by JGM
%-------------------------------------------------------------------------%


% you set these
speedThr = 60;
distThr = 40; % cm


% fixed params
T = params.dynamics.T;
Ncases = params.Ncases;
xmin = params.xmin;
xmax = params.xmax;
cntr = params.dynamics.muX0;

% get position, velocity, distance from the center, and speed
x = squeeze(S(:,:,strcmp(params.mods,'Hand-Position'),:));
v = longdata(S(:,:,strcmp(params.mods,'Hand-Velocity'),:));
cntrDist = squeeze(sqrt(sum(bsxfun(@minus,x,shiftdim(cntr,-1)).^2,2)));
mvmtSpeed = shortdata(params.Ncases,2,sqrt((v(:,1).^2 + v(:,2).^2)));


% loop through trajectories, gathering indices of center-out reaches
% for iTraj = 1:Ncases   
%     [startinds{iTraj}, endinds{iTraj}] = gatherCenterOutInds(...
%         cntrDist(iTraj,:),distThr,mvmtSpeed(iTraj,:),speedThr,T,...
%         TOPLOT,squeeze(x(iTraj,:,:)),xmin,xmax);
% end
cntrOutInds = arrayfun(@(iTraj)(gatherCenterOutInds(...
    cntrDist(iTraj,:),distThr,mvmtSpeed(iTraj,:),speedThr,T,...
    TOPLOT,squeeze(x(iTraj,:,:)),xmin,xmax)),1:Ncases,...
    'UniformOutput',false);
%%% actually, arrayfun is slightly slower than a loop here....


end