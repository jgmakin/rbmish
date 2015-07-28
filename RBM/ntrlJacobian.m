function J = ntrlJacobian(s,iMod,params)

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -now expects s to be a matrix (Ndims x Nmods) rather than a vector
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%



% init
Nmods = params.Nmods;
Ndims = params.Ndims;
sL = s(:,iMod);

% find the right Jacobian
switch params.mods{iMod}
    case {'Hand-Position','Gaze-Angle'}
        if strcmp(params.NS,'Hand-Position') || strcmp(params.NS,'Gaze-Angle')
            J = eye(Ndims);
        elseif strcmp(params.NS,'Joint-Angle')
            if Nmods == 3
                sL = s(:,1)-s(:,3);
            end
            %%%%%%%%%%%%%%%%%
            % this is a hack
            %%%%%%%%%%%%%%%%%
            J = squeeze(IKjacobianFast(sL',params));
        else
            error('unexpected params.NS');
        end
    case 'Joint-Angle'
        if strcmp(params.NS,'Hand-Position') || strcmp(params.NS,'Gaze-Angle')
            J = squeeze(FKjacobianFast(sL',params));
        elseif strcmp(params.NS,'Joint-Angle');
            J = eye(Ndims);
        else
            error('unexpected params.NS');
        end
    otherwise
        error('unrecognized modality!! -- jgm');
end



end