function posteriorCumulants = KFposteriorization(pSENSORY,Q,LDSparams,params)
% pKF = KF4PPC(unisensCmlnts,Q,LDSparams,params)
%   Takes in the sufficient-statistics structure T0 from a PPC, and filters
%   it.  This works because, at least for GTPNs, the likelihood of S given
%   the population is "Gaussian" in S---which is the requirement for the KF
%   equations to be optimal.


%-------------------------------------------------------------------------%
% Revised: 01/09/17
%   -part of Grand Revision
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

%%%% TO DO:
% (0) The idea is to parallel the (static) gaussPosteriorization
% (1) Find the functions that call this and change their behavior
% (2) functionize


% Ns
T           = Q.T;
Nexamples   = size(pSENSORY.Xpct,1);
Ndims       = params.Ndims;
Nmods       = length(params.mods);
Ntraj       = floor(Nexamples/T);


% [Nex x Ndims x Nmods] -> [Nex x Ndims*Nmods] -> [Nobsvs x T x Ntraj]
Y = reshape(pSENSORY.Xpct,[Nexamples,Ndims*Nmods]);
Y = permute(shortdata(Ntraj,3,Y),[2,3,1]);

% reshape the information matrix
if isfield(LDSparams,'SigmaYX')
    InfoYX = repmat(inv(LDSparams.SigmaYX),[1,1,T,Ntraj]);
elseif isfield(LDSparams,'InfoYX')
    InfoYX = repmat(LDSparams.InfoYX,[1,1,T,Ntraj]);
else
    % Assume no correlations b/n emissions from different populations
    % [Nex x Ndims x Ndims x Nmods] -> [Nex x Ndims*Nmods x Ndims*Nmods]
    InfoYX = zeros([Nexamples,Ndims*Nmods,Ndims*Nmods],'like',Y);
    for iMod = 1:Nmods
        inds = 1+(iMod-1)*Ndims:(iMod*Ndims);
        InfoYX(:,inds,inds) = pSENSORY.Info(:,:,:,iMod);
    end
    
    % -> [Nobsvs x Nobsvs x T x Ntraj]
    InfoYX  = permute(shortdata(Ntraj,4,InfoYX),[2,3,4,1]);
end


% Kalman filter
if isfield(Q,'U')
    U = permute(Q.U,[2,3,1]);            % put Ncases last
    KF = arrayfun(@(iCase)(KalmanFilter(setfield(LDSparams,'InfoYX',...
        InfoYX(:,:,:,iCase)),Y(:,:,iCase),'controls',U(:,:,iCase))),...
        1:size(Y,3));
else
    % KF = arrayfun(@(iCase)(KalmanFilter(setfield(LDSparams,'InfoYX',...
    %   InfoYX(:,:,:,iCase)),Y(:,:,iCase))),1:size(Y,3));
    
    XNtrpy = 0;
    fprintf('\nFiltering (Kalman)')
    for iCase = 1:Ntraj
        KF = KalmanFilter(setfield(LDSparams,'InfoYX',...
            InfoYX(:,:,:,iCase)),Y(:,:,iCase));
        fprintf('.');
        
        % store the denoised inputs
        Y(:,:,iCase) = LDSparams.C*KF.XHATMU;
        
        % you prefer to store information matrices
        CvrnYX = covX2covMX(KF.CVRNMU,LDSparams.C);
        InfoYXcell = arrayfun(@(ii)(inv(CvrnYX(:,:,ii))),1:size(CvrnYX,3),...
            'UniformOutput',false);
        InfoYX(:,:,:,iCase)  = cat(3,InfoYXcell{:});
        
        % accumulate cross entropy
        XNtrpy = XNtrpy + KF.XNtrpY;
    end
end

% reshape back to original sizes
Y = longdata(permute(Y,[3,1,2]));
if isfield(LDSparams,'muY'), Y = Y + LDSparams.muY'; end
posteriorCumulants.Xpct = reshape(Y,[Nexamples,Ndims,Nmods]);

InfoYX = longdata(permute(InfoYX,[4,1,2,3]));
posteriorCumulants.Info = zeros([Nexamples,Ndims,Ndims,Nmods],'like',Y);
for iMod = 1:Nmods
    inds = 1+(iMod-1)*Ndims:(iMod*Ndims);
    posteriorCumulants.Info(:,:,:,iMod) = InfoYX(:,inds,inds);
end

% announce
XNtrpy = XNtrpy/Ntraj;
fprintf('filtered, w/average cross entropy = %f\n\n',XNtrpy);


end
