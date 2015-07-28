function pKF = KF4PPC(LDSdata,LDSparams,name)
% pKF = KF4PPC(LDSdata,LDSparams,name)
%   Takes in the sufficient-statistics structure T0 from a PPC, and filters
%   it.  This works because, at least for GTPNs, the likelihood of S given
%   the population is "Gaussian" in S---which is the requirement for the KF
%   equations to be optimal.


%-------------------------------------------------------------------------%
% Revised: 08/27/14
%   -removed material for "noisily observed controls," since these can be
%   treated more simply as states
%   -added if statement for (fully observed) controls
% Revised: 04/30/14
%   -takes LDSdata rather than pSENSORY (frees this function from having to
%   run EFHfilter first)
% Revised: 03/10/14
%   -made some adjustments for 1D data
% Revised: 01/14/14
%   -outputs are now written directly into the relevant matrix within the
%   loop, rather than converted (via C and H, and long/shortdata) outside
% Revised: 01/06/14
%   -changed last input to "restarts" from varargin (for filterdata), and
%   -changed KFresetting appropriately
% Revised: 11/05/13
%   -now wants pSENSORY, the "posterior distribution" over X and C*X,
%   given "all the data"---i.e., just given r_t.
% Revised: 07/03/13
%   -rewrote anonymous functions to use a structure argument
%   -added case for 'resetting'
% Revised: 07/02/13
%   -drastically revised to use your KalmanFilter.m
%   -made the time-update adjustments into anonymous functions
% Revised: 07/01/13
%   -now wants T0, the sufficient stats from the (unupdated) input PPC
% Revised: 05/14/13
%   -the fxn now expects R to be just the sensory population!
% Revised: 05/11/13
%   -put prior params into params.dynamics; extracted here
% Revised: 05/09/13
%   -accommodated the wrapping of trajectories
% Revised: 05/08/13
%   -did some things
% Revised: 05/07/13
%   by JGM
% Cribbed: 05/07/13
%   from compute_kal3.m
%   by JGM & BKD
%-------------------------------------------------------------------------%


% Ns
[Ncases,Nobsvs,T] = size(LDSdata.Y);

% loop through trajectories 
C = LDSparams.C;
Y = shiftdim(LDSdata.Y,1);                          % put Ncases last
try SigmaY = shiftdim(LDSdata.SigmaY,1); 
catch ME, SigmaY = LDSparams.SigmaY; end

% Kalman filter
if isfield(LDSdata,'U')
    U = shiftdim(LDSdata.U,1);                      % put Ncases last
    KF = arrayfun(@(iCase)(KalmanFilter(setfield(LDSparams,'SigmaY',...
    SigmaY(:,:,:,iCase)),Y(:,:,iCase),U(:,:,iCase))),1:size(Y,3));
else
    % KF = arrayfun(@(iCase)(KalmanFilter(setfield(LDSparams,'SigmaY',...
    %   SigmaY(:,:,:,iCase)),Y(:,:,iCase))),1:size(Y,3));
    fprintf('\nFiltering (Kalman)')
    for iCase = 1:Ncases
        KF(iCase) = KalmanFilter(setfield(LDSparams,'SigmaY',...
            SigmaY(:,:,:,iCase)),Y(:,:,iCase));
        fprintf('.');
    end
end

% denoised output
Xpct(1:Ncases,1:Nobsvs,1:T) = shiftdim(reshape(...
    C*cat(2,KF(:).XHATMU),[Nobsvs,T,Ncases]),2);
if isfield(LDSparams,'muY')
    Xpct = Xpct + repmat(LDSparams.muY',[Ncases,1,T]);
end
Cvrn(1:Ncases,1:Nobsvs,1:Nobsvs,1:T) = shiftdim(reshape(...
    covX2covMX(cat(3,KF(:).CVRNMU),C),[Nobsvs,Nobsvs,T,Ncases]),3);
LLy = -sum([KF(:).XNtrpY])/Ncases;


% announce
fprintf('filtered, w/average LLy = %f\n\n',LLy);

% collect the outputs
pKF.Xpct = Xpct;
pKF.Cvrn = Cvrn;
pKF.params = LDSparams;
pKF.name = name;
%%%
% pKF.XHATMU = permute(cat(3,KF(:).XHATMU),[3,1,2]);
% pKF.CVRNMU = permute(cat(4,KF(:).CVRNMU),[4,1,2,3]);
%%%
%%%mat2cell(repmat(name,Nmods,1),ones(1,Nmods),length(name))';
  


end
