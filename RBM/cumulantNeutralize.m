function cumulants = cumulantNeutralize(xpctL,InfoL,params)
% cumulantNeutralize    Convert cumulants into "neutral space"
%
% USAGE:
%   posteriorCumulants = cumulantNeutralize(shatL,InfoL,params)
%
% The inputs 
%
%   xpctL (Nexamples x Ndims x Nmods) 
%
% and
%
%   InfoL  (Nexamples x Ndims x Ndims x Nmods) 
%
% are tensors containing the first two cumulants of a distribution.  But
% not all "mods" (third and fourth degrees, resp.) are in the same space.
% This function converts all the cumulants to the "neutral space" given by
% params.NS, storing the results as fields in posteriorCumulants.
%
% Note that for coordinate transformations (1Daddition, 2Daddition), this
% take in three modalities but outputs two, since two are combined.
%
% NB: shatL and InfoL needs to be in longdata format!!

%-------------------------------------------------------------------------%
% Revised: 07/06/14
%   -added code for Gaussian prior
% Revised: 07/05/14
%   -added case '1Daddition'
% Created: 07/04/14
%   by JGM
%-------------------------------------------------------------------------%


% convert non-neutral to neutral
switch params.MODEL
    case {'2Dinteg','1Dinteg','2DintegDecoupled'}
        
        % extract neutral-space and non-neutral space modalities' estimates
        sN = xpctL(:,:,strcmp(params.mods,params.NS));
        nonNSmod = params.mods{~strcmp(params.mods,params.NS)};
        snonN = xpctL(:,:,strcmp(params.mods,nonNSmod));
        
        % store
        xpctN(:,:,1) = sN;
        InfoN(:,:,:,1) = InfoL(:,:,:,strcmp(params.mods,params.NS));
        infononN = InfoL(:,:,:,strcmp(params.mods,nonNSmod));

        
        % switch on the basis of the neutral-space modality
        switch params.NS
            case 'Hand-Position'
                xpctN(:,:,2) = FK2link(snonN,params,0);
                Jinv = IKjacobianFast(sN,params);
            case 'Joint-Angle' 
                xpctN(:,:,2) = IK2link(snonN,params,0);
                Jinv = FKjacobianFast(sN,params);
            otherwise
                error('incorrect params.NS for this model');
        end
        infoJ = tensorOp(permute(infononN,[2,1,3]),permute(Jinv,[2,1,3]));
        JtrinfoJ = tensorOp(permute(Jinv,[3,1,2]),infoJ);
        InfoN(:,:,:,2) = permute(JtrinfoJ,[2,1,3]);

        % label
        srcs{1} = params.NS;
        srcs{2} = nonNSmod;
        
        
    case {'1Daddition'}
        
        % estimators and information "matrices"
        xhat = xpctL(:,:,strcmp(params.mods,'Hand-Position'));
        xInfo = InfoL(:,:,:,strcmp(params.mods,'Hand-Position'));
        thhat = xpctL(:,:,strcmp(params.mods,'Joint-Angle'));
        thInfo = InfoL(:,:,:,strcmp(params.mods,'Joint-Angle'));
        ehat = xpctL(:,:,strcmp(params.mods,'Gaze-Angle'));
        eInfo = InfoL(:,:,:,strcmp(params.mods,'Gaze-Angle'));
        
        
        % switch on the basis of the neutral-space modality
        srcs{1} = params.NS;
        switch params.NS
            case 'Hand-Position'
                xpctN(:,:,1) = xhat;
                xpctN(:,:,2) = FK2link(thhat,params,0) + ehat;
                
                InfoN(:,:,:,1) = xInfo;
                Jinv = IKjacobianFast(xhat-ehat,params);
                InfoN(:,:,:,2) =...
                    1./(1./(Jinv.*squeeze(thInfo).*Jinv) + 1./squeeze(eInfo));
                srcs{2} = 'Joint-Gaze';
                
            case 'Joint-Angle' 
                xpctN(:,:,1) = thhat;
                xpctN(:,:,2) = IK2link(xhat-ehat,params,0);
                
                InfoN(:,:,:,1) = thInfo;
                Jinv = FKjacobianFast(thhat,params);
                InfoN(:,:,:,2) =...
                    Jinv./(1./squeeze(eInfo) + 1./squeeze(xInfo)).*Jinv;
                srcs{2} = 'Hand-Gaze';
                
            case 'Gaze-Angle'
                xpctN(:,:,1) = ehat;
                xpctN(:,:,2) = xhat - FK2link(thhat,params,0);
                
                InfoN(:,:,:,1) = eInfo;
                Jinv = IKjacobianFast(xhat-ehat,params);
                InfoN(:,:,:,2) =...
                    1./(1./(Jinv.*squeeze(thInfo).*Jinv) + 1./squeeze(xInfo));
                srcs{2} = 'Joint-Hand';
                
            otherwise
                error('incorrect params.NS for this model');
        end

        
    case 'HierL2'
        
        % extract neutral-space and non-neutral space modalities' estimates
        sN = xpctL(:,:,strcmp(params.mods,params.NS));
                
        % store
        xpctN(:,:,1) = sN;
        InfoN(:,:,:,1) = InfoL(:,:,:,strcmp(params.mods,params.NS));
        srcs{1} = params.NS;
        
        % switch on the basis of the neutral-space modality
        switch params.NS
            case 'Hand-Position'
                xpctN(:,:,2) = FK2link(xpctL(:,:,strcmp(params.mods,...
                    'Joint-Angle')),params,0);
                xpctN(:,:,3) = FK2link(xpctL(:,:,strcmp(params.mods,...
                    'Joint-Angle-Left')),params,0);
                
                srcs{2} = 'Joint-Angle';
                Jinv2 = IKjacobianFast(sN,params);
                
                srcs{3} = 'Joint-Angle-Left';
                Jinv3 = Jinv2;
                    
                
            case 'Joint-Angle' 
                xpctN(:,:,2) = IK2link(xpctL(:,:,strcmp(params.mods,...
                    'Hand-Position')),params,0);
                xpctN(:,:,3) = xpctL(:,:,strcmp(params.mods,...
                    'Joint-Angle-Left'));
                
                srcs{2} = 'Hand-Position';
                Jinv2 = FKjacobianFast(sN,params);
                
                srcs{3} = 'Joint-Angle-Left';
                I = eye(size(sN,2)); Ivec = I(:)';
                Jinv3 = reshape(Ivec(ones(size(sN,1),1),:),...
                    [size(sN,1),size(sN,2),size(sN,2)]);
                
                
            case 'Joint-Angle-Left'
                xpctN(:,:,2) = IK2link(xpctL(:,:,strcmp(params.mods,...
                    'Hand-Position')),params,0);
                xpctN(:,:,3) = xpctL(:,:,strcmp(params.mods,...
                    'Joint-Angle'));
                
                srcs{2} = 'Hand-Position';
                Jinv2 = FKjacobianFast(sN,params);                
                
                srcs{3} = 'Joint-Angle';
                I = eye(size(sN,2)); Ivec = I(:)';
                Jinv3 = reshape(Ivec(ones(size(sN,1),1),:),...
                    [size(sN,1),size(sN,2),size(sN,2)]);
                
                
            otherwise
                error('incorrect params.NS for this model');
        end
        info = InfoL(:,:,:,strcmp(params.mods,srcs{2}));
        infoJ = tensorOp(permute(info,[2,1,3]),permute(Jinv2,[2,1,3]));
        JtrinfoJ = tensorOp(permute(Jinv2,[3,1,2]),infoJ);
        InfoN(:,:,:,2) = permute(JtrinfoJ,[2,1,3]);
        
        info = InfoL(:,:,:,strcmp(params.mods,srcs{3}));
        infoJ = tensorOp(permute(info,[2,1,3]),permute(Jinv3,[2,1,3]));
        JtrinfoJ = tensorOp(permute(Jinv3,[3,1,2]),infoJ);
        InfoN(:,:,:,3) = permute(JtrinfoJ,[2,1,3]);
        
        
    case {'2Dtwoarms'}
        
        % everything is "neutral"
        xpctN = xpctL;
        InfoN = InfoL;
        
        % label
        srcs{1} = params.NS;
        srcs{2} = params.mods{~strcmp(params.mods,params.NS)};
        
    otherwise
        error('unknown model for cumulantNeutralize -- jgm\n')
end

% add prior distribution, if there is one
if isfield(params,'p0')
    xpct0 = params.p0.mu;
    Info0 = inv(params.p0.cov);
    
    xpctN(:,:,end+1) = xpct0(:,ones(size(xpctN,1),1))';
    InfoN(:,:,:,end+1) = permute(Info0(:,:,ones(size(xpctN,1),1)),[3,1,2]);
    srcs{end+1} = 'prior';
end


% store
cumulants.Xpct = xpctN;
cumulants.Info = InfoN;
cumulants.srcs = srcs;


end
